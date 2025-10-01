# LikelihoodProfiler
*A package for practical identifiability analysis and confidence intervals (CI) estimation using the profile likelihood approach. The package provides a unified interface for various profile likelihood methods, including optimization-based and integration-based profiles, CI endpoints search, and more.*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://insysbio.github.io/LikelihoodProfiler.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://insysbio.github.io/LikelihoodProfiler.jl/latest)
[![Build Status](https://github.com/insysbio/LikelihoodProfiler.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/insysbio/LikelihoodProfiler.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/insysbio/LikelihoodProfiler.jl/graph/badge.svg)](https://codecov.io/gh/insysbio/LikelihoodProfiler.jl)

## Installation

In Julia terminal run the following command:

```julia
import Pkg; Pkg.add("LikelihoodProfiler")
```

## Breaking changes and new features

**NOTE:** version 1.2.0 introduces profiling **functions of parameters** while keeping a single `ProfileLikelihoodProblem` interface. To make this possible (and more consistent with `Optimization.jl`), we made a few **breaking adjustments**.

  • `idxs` moved from `solve(…; idxs=…)` to the problem constructor. You now select which parameters to profile when building the `ProfileLikelihoodProblem`.

  • `profile_range` (positional) was argument replaced by `profile_lower / profile_upper` (keywords). This mirrors `Optimization.jl`’s style and avoids confusion with `optprob.lb/optprob.ub`.

### Migrating your code

1) Parameter profiling

Before (≤ 1.1.x):
```julia
plprob = ProfileLikelihoodProblem(optprob, optpars, profile_range)  # e.g. [(-5,5), (-2,2), …]
sol = solve(plprob, method; idxs=[1,3])
```

After (≥ 1.2.0):
```julia
# Pass indices and profile bounds at construction.
plprob = ProfileLikelihoodProblem(optprob, optpars;
                                  idxs=[1,3],
                                  profile_lower=[-5.0, -1.0],
                                  profile_upper=[ 2.0,  4.0])

sol = solve(plprob, method)  # no `idxs` here anymore
```
  • If you omit `profile_lower/profile_upper`, and your `OptimizationProblem` has `lb/ub`, they will be used (sliced by `idxs`).

1) Function profiling (new)
```julia
g1 = OptimizationFunction((θ,p)->θ[1] + θ[2])
g2 = OptimizationFunction((θ,p)->θ[2] - θ[3])

plprob = ProfileLikelihoodProblem(optprob, optpars, [g1,g2];
                                  profile_lower=-2.0,    # scalars expand per target
                                  profile_upper= 2.0)

sol = solve(plprob, method)
```

## Getting started with LikelihoodProfiler

To define a profile likelihood problem `ProfileLikelihoodProblem` in LikelihoodProfiler, you should provide the objective function (usually negative log likelihood) and the optimal values of the parameters that correspond to the minimum of the objective function. LikelihoodProfiler relies on the `Optimization.jl` interface, and `ProfileLikelihoodProblem` is built on top of the `OptimizationProblem` defined in `Optimization.jl`. This can be best illustrated by an example.

First we define the `OptimizationProblem` and solve it with the preferred optimizer to obtain the optimal values of the parameters. 

```julia
using Optimization, ForwardDiff

# objective function
rosenbrock(x,p) = (1.0 - x[1])^2 + 100.0*(x[2] - x[1]^2)^2

# initial values
x0 = zeros(2)

# solving optimization problem
optf = OptimizationFunction(rosenbrock, AutoForwardDiff())
optprob = OptimizationProblem(optf, x0)
sol = solve(optprob, Optimization.LBFGS())
```

### Profile likelihood problem interface

To define the `ProfileLikelihoodProblem`, we need the `OptimizationProblem` and the optimal values of the parameters. We can also set the profiling domain with the `profile_lower`, `profile_upper` arguments, indicies of parameters to profile with `idxs` and the `threshold`, which is the confidence level required to estimate confidence intervals. Please consult `?ProfileLikelihoodProblem` on the details of the interface.

```julia
using LikelihoodProfiler, Plots

# optimal values of the parameters
optpars = sol.u

# profile likelihood problem
plprob = ProfileLikelihoodProblem(optprob, optpars; profile_lower = -10., profile_upper=10., threshold = 4.0)
```

### Profile likelihood methods

LikelihoodProfiler provides a range of methods to profile likelihood functions and explore practical identifiability. The most common and simple "profiler" is the `OptimizationProfiler` method. It is based on stepwise re-optimization of the likelihood function with the constraint on the parameter (or function) of interest. We define the method and run the `solve` procedure. Please consult `?solve` on the details of the interface.

```julia
method = OptimizationProfiler(optimizer = Optimization.LBFGS(), stepper = FixedStep(; initial_step=0.15))
sol = solve(plprob, method)
plot(sol, size=(800,300), margins=5Plots.mm)
```
![Rosenbrock optimization-based profile](https://github.com/insysbio/LikelihoodProfiler.jl/blob/master/docs/assets/rosenbrock_optimization.png)

The same `solve` interface can be used with other profiling methods. For example, a more advanced way to compute profiles is proposed by `IntegrationProfiler`. It obtains the profiles as solutions to the differential equation system. To solve this internally generated system, we need to provide a differential equations solver (`integrator`). 

```julia
using OrdinaryDiffEq

method = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.3,), matrix_type = :hessian)
sol = solve(plprob, method)
plot(sol, size=(800,300), margins=5Plots.mm)
```
![Rosenbrock integration-based profile](https://github.com/insysbio/LikelihoodProfiler.jl/blob/master/docs/assets/rosenbrock_integration.png)

Likelihood profiling is mostly performed to assess if the profile has intersections with the given confidence level, hence if the parameter (or function of parameters) has finite confidence interval. Another approach to the problem of practical identifiability is to compute these intersections (endpoints of the confidence interval (CI)) without restoring the full shape of the profile. One of such methods is implemented in `CICOProfiler`. It estimates CI endpoints with an optimization procedure without following the exact trajectory of the profile. 

```julia
using CICOBase

method = CICOProfiler(optimizer = :LN_NELDERMEAD, scan_tol = 1e-4)
sol = solve(plprob, method)
plot(sol, size=(800,300), margins=5Plots.mm)
```
![Rosenbrock CICO profile](https://github.com/insysbio/LikelihoodProfiler.jl/blob/master/docs/assets/rosenbrock_cico.png)

## Related packages

Other implementations of the profile likelihood approach in Julia include:
- [ProfileLikelihood.jl](https://github.com/DanielVandH/ProfileLikelihood.jl) implements fixed-step optimization-based profiles and supports bivariate profile likelihoods.
- [InformationalGeometry.jl](https://github.com/RafaelArutjunjan/InformationGeometry.jl) implements various methods to study likelihood functions (including profile likelihood) using the tools of differential geometry.

There are also well-known profile likelihood implementations in other languages, namely: [Data2Dynamics](https://github.com/Data2Dynamics/d2d), [dMod](https://github.com/dkaschek/dMod/), [pyPESTO](https://github.com/ICB-DCM/pyPESTO), [sbioparametersci](https://www.mathworks.com/help/simbio/ref/sbioparameterci.html)

## Citation

Borisov I., Metelkin E. An Algorithm for Practical Identifiability Analysis and Confidence Intervals Evaluation Based on Constrained Optimization. 2018. October. ICSB2018. https://doi.org/10.13140/RG.2.2.18935.06563
