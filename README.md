# LikelihoodProfiler
*A package for practical identifiability analysis and confidence intervals (CI) estimation using profile likelihood approach. The package provides a single interface to different profile likelihood methods, including optimization-based, integration-based profiles, CI endpoints search and more.*

#### !!! LikelihoodProfiler package underwent a huge redevelopment and now provides a single interface to different profile likelihood methods. If you are looking for a the CICO algorithm (initially implemented in LikelihoodProfiler package) you can use it through the new LikelihoodProfiler interface or directly with [CICOBase.jl package](https://github.com/insysbio/CICOBase.jl).

#### !!! The new interface is in the process of development. Feedback and contributions are very welcome!  

## Installation

In Julia terminal run the following command:

```julia
import Pkg; Pkg.add("LikelihoodProfiler")
```

## Getting started with LikelihoodProfiler

To define a profile likelihood problem `PLProblem` in LikelihoodProfiler you should provide the objective function (usually negative log likelihood) and the optimal values of the parameters which correspond to the minimum of the objective function. LikelihoodProfiler relies on `Optimization.jl` interface and `PLProblem` is build on top of the `OptimizationProblem` defined in `Optimization.jl`. This can be best illustrated by an example.

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

To define the `PLProblem` we need the `OptimizationProblem` and the optimal values of the parameters. Also we can set profiling domain with `profile_range` argument and `threshold` which is the confidence level required to estimate confidence intervals. Please, consult `?PLProblem` on the details of the interface.

```julia
using LikelihoodProfiler, Plots

# optimal values of the parameters
optpars = sol.u

# profile likelihood problem
plprob = PLProblem(optprob, optpars, (-10.,10.); threshold = 4.0)
```

### Profile likelihood methods

LikelihoodProfiler provides a range of methods to profile likelihood functions and explore practical identifiability. The most common and simple "profiler" is the `OptimizationProfiler` method. It is based on the stepwise re-optimization of the likelihood function with the constraint on the parameter (or function) of interest. We define the method and run `profile` procedure. Please, consult `?profile` on the details of the interface.

```julia
method = OptimizationProfiler(optimizer = Optimization.LBFGS(), stepper = FixedStep(; initial_step=0.15))
sol = profile(plprob, method)
plot(sol, size=(800,300), margins=5Plots.mm)
```
![Rosenbrock optimization-based profile](https://github.com/insysbio/LikelihoodProfiler.jl/blob/dev/docs/assets/rosenbrock_optimization.png)

The same `profile` interface can be used with other profiling methods. For example a more advanced way to compute profiles is proposed by `IntegrationProfiler`. It obtains the profiles as the solution to the differential equation system. In order to solve this internally generated system we need to provide a differential equations solver (`integrator`). 

```julia
using OrdinaryDiffEq

method = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.3,), matrix_type = :hessian)
sol = profile(plprob, method)
plot(sol, size=(800,300), margins=5Plots.mm)
```
![Rosenbrock integration-based profile](https://github.com/insysbio/LikelihoodProfiler.jl/blob/dev/docs/assets/rosenbrock_integration.png)

Likelihood profiling is mostly performed to assess if the profile has intersections with the given confidence level, hence if the parameter (or function of parameters) has finite confidence interval. Another approach to the problem of practical identifiability is to compute these intersections (endpoints of the confidence interval (CI)) without restoring the full shape of the profile. One of such methods is implemented in `CICOProfiler`. It estimates CI endpoints with an optimization procedure without following the exact trajectory of the profile. 

```julia
using CICOBase

method = CICOProfiler(optimizer = :LN_NELDERMEAD, scan_tol = 1e-4)
sol = profile(plprob, method)
plot(sol, size=(800,300), margins=5Plots.mm)
```
![Rosenbrock CICO profile](https://github.com/insysbio/LikelihoodProfiler.jl/blob/dev/docs/assets/rosenbrock_cico.png)

## Related packages

Other implementations of the profile likelihood approach in Julia include:
- [ProfileLikelihood.jl](https://github.com/DanielVandH/ProfileLikelihood.jl) implements fixed-step optimization-based profiles. Also supports bivariate profile likelihoods.
- [InformationalGeometry.jl](https://github.com/RafaelArutjunjan/InformationGeometry.jl) implements various methods to study likelihood functions (including profile likelihood) with the tools of differential geometry.

Also there are well-known profile likelihood implementations in other languages, namely: [Data2Dynamics](https://github.com/Data2Dynamics/d2d), [dMod](https://github.com/dkaschek/dMod/), [pyPESTO](https://github.com/ICB-DCM/pyPESTO), [sbioparametersci](https://www.mathworks.com/help/simbio/ref/sbioparameterci.html)

## Citation

Borisov I., Metelkin E. An Algorithm for Practical Identifiability Analysis and Confidence Intervals Evaluation Based on Constrained Optimization. 2018. October. ICSB2018. https://doi.org/10.13140/RG.2.2.18935.06563
