# 1. Import
import logging
import tempfile
from pprint import pprint
from collections.abc import Sequence
from typing import Union

import amici
import matplotlib as mpl
import numpy as np
import petab

import pypesto.optimize as optimize
import pypesto.petab
import pypesto.profile as profile
import pypesto.sample as sample
import pypesto.store as store
import pypesto.visualize as visualize
import pypesto.visualize.model_fit as model_fit

from pypesto.profile import calculate_approximate_ci, chi2_quantile_to_ratio
from pypesto.result import Result

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np

# Set seed for reproducibility
np.random.seed(1912)

# 2. Import the model
model_name = "Boehm_JProteomeRes2014"
petab_yaml = f"../Benchmark-Models-PEtab/Benchmark-Models/{model_name}/{model_name}.yaml"

petab_problem = petab.Problem.from_yaml(petab_yaml)
importer = pypesto.petab.PetabImporter(petab_problem)
problem = importer.create_problem(verbose=False)

# 3. Setup optimizers
optimizer_options = {"maxiter": 1e4, "fatol": 1e-12, "frtol": 1e-12}
optimizer = optimize.FidesOptimizer(options=optimizer_options, verbose=logging.WARN)
optimizer_scipy_lbfgsb = optimize.ScipyOptimizer(method="L-BFGS-B")
optimizer_scipy_powell = optimize.ScipyOptimizer(method="Powell")

opt_options = optimize.OptimizeOptions()

n_starts = 20  # usually a value >= 100 should be used
engine = pypesto.engine.SingleCoreEngine()

# 4. Minimize
result = optimize.minimize(
    problem=problem,
    optimizer=optimizer_scipy_lbfgsb,
    n_starts=n_starts,
    engine=engine,
    options=opt_options,
)

# 5. Profile
result = profile.parameter_profile(
    problem=problem,
    result=result,
    optimizer=optimizer_scipy_lbfgsb,
    engine=engine,
)

# 6. Print the CIs
def profile_cis(
    result: Result,
    confidence_level: float = 0.95,
    profile_indices: Sequence[int] = None,
    profile_list: int = 0,
    color: Union[str, tuple] = "C0",
    show_bounds: bool = False,
    ax: matplotlib.axes.Axes = None,
) -> matplotlib.axes.Axes:
    """
    Plot approximate confidence intervals based on profiles.

    Parameters
    ----------
    result:
        The result object after profiling.
    confidence_level:
        The confidence level in (0,1), which is translated to an approximate
        threshold assuming a chi2 distribution, using
        `pypesto.profile.chi2_quantile_to_ratio`.
    profile_indices:
        List of integer values specifying which profiles should be plotted.
        Defaults to the indices for which profiles were generated in profile
        list `profile_list`.
    profile_list:
        Index of the profile list to be used.
    color:
        Main plot color.
    show_bounds:
        Whether to show, and extend the plot to, the lower and upper bounds.
    ax:
        Axes object to use. Default: Create a new one.
    """
    # extract problem
    problem = result.problem
    # extract profile list
    profile_list = result.profile_result.list[profile_list]

    if profile_indices is None:
        profile_indices = [ix for ix, res in enumerate(profile_list) if res]

    if ax is None:
        fig, ax = plt.subplots()

    confidence_ratio = chi2_quantile_to_ratio(confidence_level)

    # calculate intervals
    res = []
    intervals = []
    for i_par in range(problem.dim_full):
        if i_par not in profile_indices:
            continue
        xs = profile_list[i_par].x_path[i_par]
        ratios = profile_list[i_par].ratio_path
        lb, ub = calculate_approximate_ci(
            xs=xs, ratios=ratios, confidence_ratio=confidence_ratio
        )
        res.append((problem.x_names[i_par], (lb, ub)))
        intervals.append((lb, ub))

    x_names = [problem.x_names[ix] for ix in profile_indices]

    for ix, (lb, ub) in enumerate(intervals):
        ax.plot([lb, ub], [ix + 1, ix + 1], marker="|", color=color)

    parameters_ind = np.arange(1, len(intervals) + 1)
    ax.set_yticks(parameters_ind)
    ax.set_yticklabels(x_names)
    ax.set_ylabel("Parameter")
    ax.set_xlabel("Parameter value")

    if show_bounds:
        lb = problem.lb_full[profile_indices]
        ax.plot(lb, parameters_ind, "k--", marker="+")
        ub = problem.ub_full[profile_indices]
        ax.plot(ub, parameters_ind, "k--", marker="+")

    return res, fig, ax

res, fig, ax = profile_cis(result)

for (param, (lb, ub)) in res:
    print(param, ", CI is ", (lb, ub))

fig.savefig("boehm_pypesto.svg")

