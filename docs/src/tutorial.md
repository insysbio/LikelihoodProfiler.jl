## Getting started with LikelihoodProfiler

To define a profile likelihood problem `PLProblem` in LikelihoodProfiler you should provide the objective function (usually negative log likelihood) and the optimal values of the parameters which correspond to the minimum of the objective function. LikelihoodProfiler relies on `Optimization.jl` interface and `PLProblem` is build on top of the `OptimizationProblem` defined in `Optimization.jl`. This can be best illustrated by an example.

First we define the `OptimizationProblem` and solve it with the preferred optimizer to obtain the optimal values of the parameters. 

```@example example-1
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

```@example example-1
using LikelihoodProfiler, Plots

# optimal values of the parameters
optpars = sol.u

# profile likelihood problem
plprob = PLProblem(optprob, optpars, (-10.,10.); threshold = 4.0)
```

### Profile likelihood methods

LikelihoodProfiler provides a range of methods to profile likelihood functions and explore practical identifiability. The most common and simple "profiler" is the `OptimizationProfiler` method. It is based on the stepwise re-optimization of the likelihood function with the constraint on the parameter (or function) of interest. We define the method and run `profile` procedure. Please, consult `?profile` on the details of the interface.

```@example example-1
method = OptimizationProfiler(optimizer = Optimization.LBFGS(), stepper = FixedStep(; initial_step=0.15))
sol = profile(plprob, method)
plot(sol, size=(800,300), margins=5Plots.mm)
```

The same `profile` interface can be used with other profiling methods. For example a more advanced way to compute profiles is proposed by `IntegrationProfiler`. It obtains the profiles as the solution to the differential equation system. In order to solve this internally generated system we need to provide a differential equations solver (`integrator`). 

```@example example-1
using OrdinaryDiffEq

method = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.3,), matrix_type = :hessian)
sol = profile(plprob, method)
plot(sol, size=(800,300), margins=5Plots.mm)
```
