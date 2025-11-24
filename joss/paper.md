---
title: "LikelihoodProfiler.jl: Unified profile-likelihood workflows for identifiability and confidence intervals"
tags:
  - profile likelihood
  - identifiability analysis
  - parameter estimation
  - systems biology
  - quantitative systems pharmacology
  - mathematical modeling
authors:
  - name: "Ivan Borisov"
    orcid: 0000-0002-2887-1903
    affiliation: 1
  - name: "Aleksander Demin"
    affiliation: 2
  - name: "Evgeny Metelkin"
    orcid: 0000-0001-6612-5373
    affiliation: 1
affiliations:
  - name: "InSysBio CY LTD, Limassol, Cyprus"
    index: 1
  - name: "National Research University Higher School of Economics, Moscow, Russia"
    index: 2
date: 12 August 2025
bibliography: paper.bib
---

## Summary

Practical identifiability addresses how well a mechanistic model is determined by available experimental data. Profile likelihood-based methods are widely used to explore practical identifiability and often serve as a proxy for structural identifiability analysis [@Heinrich2025]. Moreover, these methods can also be generalized to states and predictions analysis. This versatility makes profile likelihood an essential component in model development.

The paper presents `LikelihoodProfiler.jl`, an open-source Julia package for conducting profile likelihood-based identifiability and predictability analysis.

## Statement of Need

Current profile likelihood implementations are often tied to specific modeling frameworks or software environments. In addition, different profile likelihood algorithms are spread across separate tools, forcing users to switch between languages and interfaces.
`LikelihoodProfiler.jl` addresses these limitations by providing:
- Unified interface to multiple profile likelihood algorithms.
- Compatibility with common modeling standards (Heta [@Metelkin2021], SBML, PEtab [@Persson2025]).
- Integration with Julia’s SciML [@Rackauckas2017] for efficient computation and extensibility.

## Features and Methodologies

`LikelihoodProfiler.jl` provides a unified interface for various profile likelihood methods, including optimization-based **OptimizationProfiler** and integration-based profiles **IntegrationProfiler**, CI endpoints search **CICOProfiler**, and more.

All methods leverage a `CommonSolve` interface [@Rackauckas2017] (`CommonSolve.solve()`), supporting global settings for parallelization and verbosity. Profiling results are directly visualizable using the `Plots.jl` package and exportable as DataFrames for further analysis.

## Demonstrative Example: JAK/STAT Signaling Pathway Model

`LikelihoodProfiler.jl`’s functionality and interfaces are demonstrated using the JAK/STAT signaling pathway model [@Boehm2014], which consists of 8 states and 9 parameters. The model and experimental data were sourced from the Benchmark-Models-PEtab repository [@petab_benchmark_collection] and imported through the `PEtab.jl` package [@Persson2025].

```julia
using LikelihoodProfiler, CICOBase, PEtab, Plots
petab_model = PEtabModel("Boehm_JProteomeRes2014.yaml")
petab_problem = PEtabODEProblem(petab_model)
```

To study the identifiability of the JAK/STAT model parameters, we first construct a `ProfileLikelihoodProblem` and select an appropriate profiling method.
A `ProfileLikelihoodProblem` is defined by providing (i) the objective function (typically the negative log-likelihood) and (ii) the initial parameter values corresponding to the optimum of this objective. `LikelihoodProfiler.jl` builds on the `Optimization.jl` interface [@vaibhav_kumar_dixit_2023_7738525], and  `ProfileLikelihoodProblem` wraps an `OptimizationProblem` defined in `Optimization.jl`. In addition, `ProfileLikelihoodProblem` allows users to specify optional arguments common across profiling methods, such as the indices of parameters to profile, upper and lower bounds, and other options.

```julia
using Optimization, OrdinaryDiffEq
optprob = OptimizationProblem(petab_problem)
param_profile_prob = ProfileLikelihoodProblem(optprob, get_x(petab_problem); idxs=1:3)
```
`LikelihoodProfiler.jl` offers a suite of methods for solving the `ProfileLikelihoodProblem`. In this example, we use the Hessian-free variant of the `IntegrationProfiler` [@Chen2002], which approximates the likelihood profile and performs re-optimization after each ODE solver step to prevent divergence from the true profile trajectory. Each profiling method offers several configurable options.
```julia
alg = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.5, reltol=1e-3, abstol=1e-4),
  matrix_type = :identity, gamma=0., reoptimize=true, 
  optimizer = Optimization.LBFGS(),optimizer_opts=(maxiters=10000,))

sol_param = solve(param_profile_prob, alg, verbose=true)
```


An alternative approach, implemented in `CICOProfiler`, estimates the confidence intervals (CI) endpoints directly — without reconstructing the full profile — by solving a constrained optimization problem [@Borisov2020]. This method is often more efficient when only the CI is required.

```julia
alg_cico = CICOProfiler(optimizer = :LN_NELDERMEAD, scan_tol = 1e-10)
```
Users can get estimated confidence intervals using the `endpoints(sol_param)` function or visualize profile curves with `plot()`.

![Parameter profiles](param_profiles.png)

![CICO profiles](cico_profiles.png)

The same algorithms can also be applied to arbitrary functions of parameters. This generalization of profile likelihood concept [@Kreutz2012] is used to study model reparametrization or perform predictability analysis. This functionality is available within the same interface: users define functions of parameters and provide them as the third argument to the problem. As an illustration, we define a function that returns the observable value at selected time points:

```julia
function pSTAT5A_rel_obs(x, p, t)
  sol = get_odesol(x, petab_problem)(t)
  specC17 = 0.107
  return (100 * sol[7] + 200 * sol[6] * specC17) / (sol[7] + sol[1] * specC17 + 2 * sol[6] * specC17)
end
pSTAT5A_rel_optf(t) = OptimizationFunction((x, p) -> pSTAT5A_rel_obs(x, p, t), AutoFiniteDiff())

times = [2.5, 60.0, 200.0]
func_profile_prob = ProfileLikelihoodProblem(optprob, get_x(petab_problem), [pSTAT5A_rel_optf(t) for t in times];
  profile_lower=0.0, profile_upper=120.0)

sol_func = solve(func_profile_prob, alg)
```
![Function profiles](func_profiles.png)

Prediction confidence interval endpoints can be estimated at a larger set of time points, and smooth prediction bands can then be plotted.

![Prediction profile](pred_profile.png)

## Implementation and Extensibility

All profiling methods benefit from the unified interface provided by `LikelihoodProfiler.jl`:
- Integration with `SciML` packages gives users access to a wide range of optimizers, differential equation solvers, and AD backends, enabling efficient profiling configurations.
- Compatibility with `Heta`, `PEtab` and `SBML` formats broadens the accessibility of the package across different modeling frameworks.
- `solve()` interface provided by `CommonSolve.jl` provides unified access various profiling methods which can be used both for parameters and for functions identifiability analysis
- A common parallelization setup, controlled via the `parallel_type` argument in the `solve()` function, is supported across all methods and can significantly accelerate computations.
- The interface facilitates integration of new profiling methods and stepping algorithms.

Future work will include adding new methods of parameters, functions and predictions profiling and enabling adaptive switching between strategies.

## Availability

LikelihoodProfiler.jl is open-source and available at: https://github.com/insysbio/LikelihoodProfiler.jl

The package is registered in Julia and can be installed from the Julia REPL using:
```julia
import Pkg; Pkg.add("LikelihoodProfiler")
```

Tutorials and documentation are available at: https://insysbio.github.io/LikelihoodProfiler.jl/stable/

## Related packages

Julia packages, such as `ProfileLikelihood.jl` and `InformationGeometry.jl`, offer related functionality but do not provide multiple profiling methods through a unified interface. Outside Julia, tools like `Data2Dynamics (MATLAB)`, `dMod (R)`, and `pyPESTO (Python)` combine parameter estimation and profile likelihoods. However, these frameworks are typically tied to their respective modeling environments and are not designed to support identifiability analysis independently of the model-development workflow.

## References
