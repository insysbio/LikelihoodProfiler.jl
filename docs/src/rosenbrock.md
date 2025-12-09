## Getting Started with LikelihoodProfiler

`LikelihoodProfiler.jl` provides tools for computing **profile likelihoods**, **confidence intervals**, and **prediction profiles** for maximum-likelihood estimation (MLE) problems.  
To define a `ProfileLikelihoodProblem`, you need two ingredients:

1. an **objective function** (typically a negative log-likelihood or a least-squares loss), and  
2. the **optimal parameter values**, i.e., the parameter vector that minimizes the objective.

Profile likelihood analysis requires evaluating the objective function under parameter perturbations. To define the objective function and its associated optimization problem `LikelihoodProfiler.jl` builds directly on the `Optimization.jl` interface. Internally, every `ProfileLikelihoodProblem` wraps an `OptimizationProblem`, so anything you can optimize with `Optimization.jl` can be profiled here.

Below we illustrate the basic workflow using a simple artificial objective function: the [Rosenbrock function](https://en.wikipedia.org/wiki/Rosenbrock_function).  
More realistic applications (ODE models, statistical models) are shown in the Tutorials section.

First we define the `OptimizationProblem` and solve it with the preferred optimizer to obtain the optimal values of the parameters. 

### Define and solve the optimization problem

We begin by defining the objective function and solving the optimization problem to obtain the optimal values of the parameters:

```@example rosenbrock-1
using LikelihoodProfiler, OptimizationLBFGSB, OrdinaryDiffEq, CICOBase
using Plots

# objective function
rosenbrock(x,p) = (1.0 - x[1])^2 + 100.0*(x[2] - x[1]^2)^2

# initial values
x0 = zeros(2)

# solving optimization problem
optf = OptimizationFunction(rosenbrock, AutoForwardDiff())
optprob = OptimizationProblem(optf, x0)
sol = solve(optprob, LBFGSB())
```

### Profile Likelihood Problem Interface

To construct a `ProfileLikelihoodProblem`, we provide the `OptimizationProblem` together with the optimal parameter values obtained from solving it. We can also set the profiling domain with the `profile_lower`, `profile_upper` arguments, indicies of parameters to profile with `idxs` and the `threshold`, which is the confidence level required to estimate confidence intervals. Please consult `?ProfileLikelihoodProblem` on the details of the interface.

```@example rosenbrock-1
# optimal values of the parameters
optpars = sol.u

# profile likelihood problem
plprob = ProfileLikelihoodProblem(optprob, optpars; profile_lower = -10., profile_upper=10., threshold = 4.0)
```

### Profile Likelihood Methods

`LikelihoodProfiler` provides a range of methods to profile likelihood functions and explore practical identifiability. All methods use the same [`solve`](@ref solve) interface but differ in how the profiles and the associated confidence intervals are computed. The `solve` interface exposes a set of common options shared across all profiling methods—such as `parallel_type` for parallel execution, `verbose` for controlling output, and `reoptimize_init` to trigger re-optimization of the initial parameter values before profiling begins.

For a more detailed description of each method, see the [Profile Likelihood Methods](@ref profile_likelihood_methods) section.

#### OptimizationProfiler

The most direct profiling method is [`OptimizationProfiler`](@ref optimization_based_profiles). It computes the profile by taking a sequence of fixed steps in the profiled parameter and performing a re-optimization of all remaining parameters at each step. This produces a discrete set of points that approximate the profile curve.

```@example rosenbrock-1
meth_opt = OptimizationProfiler(optimizer = LBFGSB(), stepper = FixedStep(; initial_step=0.15))
sol1 = solve(plprob, meth_opt)
plot(sol1, size=(800,300), margins=5Plots.mm)
```

#### IntegrationProfiler

A more advanced technique is implemented in [`IntegrationProfiler`](@ref integration_based_profiles). Instead of repeatedly re-optimizing parameters, this method internally formulates a differential equation system whose solution traces the profile likelihood trajectory. To integrate this system, the user must provide a differential equation solver (`integrator`).

```@example rosenbrock-1
meth_integ = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.3,), matrix_type = :hessian)
sol2 = solve(plprob, meth_integ)
plot(sol2, size=(800,300), margins=5Plots.mm)
```

#### CICOProfiler

Often, the primary goal of likelihood profiling is to determine whether the profile intersects the confidence threshold—i.e., whether the parameter has a finite confidence interval.
[`CICOProfiler`](@ref cico_profiles) focuses on this directly: it estimates the confidence interval endpoints without reconstructing the full profile curve.

```@example rosenbrock-1
meth_cico = CICOProfiler(optimizer = :LN_NELDERMEAD, scan_tol = 1e-4)
sol3 = solve(plprob, meth_cico)
plot(sol3, size=(800,300), margins=5Plots.mm)
```

### Profile Likelihood Solution
A `ProfileLikelihoodSolution` stores more than just the profile curves (see [Solution Interface](@ref solution_interface)). It also contains the estimated confidence-interval endpoints and identification retcodes, which indicate whether each parameter (or function of parameters) is practically identifiable.
These values can be accessed directly:

```@example rosenbrock-1
retcodes(sol3)
endpoints(sol3)
```
