# LikelihoodProfiler
*A package for [practical identifiability analysis](https://en.wikipedia.org/wiki/Identifiability_analysis) and confidence intervals estimation using profile likelihood approach. The package provides a single interface to different profile likelihood methods.*

#### !!! LikelihoodProfiler package underwent a huge redevelopment and now provides a single interface to different profile likelihood methods. If you are looking for a the CICO algorithm (initially implemented in LikelihoodProfiler package) you can use it through the new LikelihoodProfiler interface or directly with [CICO.jl package](https://github.com/insysbio/CICO.jl).

#### !!! The new interface is at an early stage. Feedback and contributions are very welcome!  

## Installation

```julia
julia> ]

(v1.11) pkg> add LikelihoodProfiler
```

## Getting started with LikelihoodProfiler

#### !!! This section is the draft of the documentation. Detailed documentation of the interface and profiling methods is coming soon.

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
using LikelihoodProfiler

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
plot(sol, size=(800,300))
```
![Rosenbrock optimization-based profile](https://github.com/insysbio/LikelihoodProfiler.jl/blob/develop/docs/assets/rosenbrock_optimization.png)

The same `profile` interface can be used with other profiling methods. For example a more advanced way to compute profiles is proposed by `IntegrationProfiler`. It obtains the profiles as the solution to the differential equation system. In order to solve this internally generated system we need to provide a differential equations solver (`integrator`). 

```julia
using OrdinaryDiffEq

method = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.3,), matrix_type = :hessian)
sol = profile(plprob, method)
plot(sol, size=(800,300))
```
![Rosenbrock integration-based profile](https://github.com/insysbio/LikelihoodProfiler.jl/blob/develop/docs/assets/rosenbrock_integration.png)

Likelihood profiling is mostly performed to assess if the profile has intersections with the given confidence level, hence if the parameter (or functions) has finite confidence intervals. Another approach to the problem of practical identifiability is to compute these intersections (endpoints of the confidence interval (CI)) without restoring the full shape of the profile. On of such methods is implemented in `CICOProfiler`. It estimates CI endpoints with an optimization procedure without following the exact trajectory of the profile. 

#### (TODO)

## License

[MIT](./LICENSE) Public license

## How to cite

**Borisov I, Metelkin E** (2020) *Confidence intervals by constrained optimization—An algorithm and software package for practical identifiability analysis in systems biology.* PLoS Comput Biol 16(12): e1008495. Ref: <https://doi.org/10.1371/journal.pcbi.1008495>
