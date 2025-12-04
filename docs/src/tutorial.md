## Getting started with LikelihoodProfiler

To define a profile likelihood problem `ProfileLikelihoodProblem` in LikelihoodProfiler, you should provide the objective function (usually negative log likelihood) and the optimal values of the parameters that correspond to the minimum of the objective function. LikelihoodProfiler relies on the `Optimization.jl` interface, and `ProfileLikelihoodProblem` is built on top of the `OptimizationProblem`. This can be best illustrated by an example.

First we define the `OptimizationProblem` and solve it with the preferred optimizer to obtain the optimal values of the parameters. 

```@example example-1
using OptimizationLBFGSB, ForwardDiff

# objective function
rosenbrock(x,p) = (1.0 - x[1])^2 + 100.0*(x[2] - x[1]^2)^2

# initial values
x0 = zeros(2)

# solving optimization problem
optf = OptimizationFunction(rosenbrock, AutoForwardDiff())
optprob = OptimizationProblem(optf, x0)
sol = solve(optprob, LBFGSB())
```

### Profile likelihood problem interface

To define the `ProfileLikelihoodProblem`, we need the `OptimizationProblem` and the optimal values of the parameters. We can also set the profiling domain with the `profile_lower`, `profile_upper` arguments, indicies of parameters to profile with `idxs` and the `threshold`, which is the confidence level required to estimate confidence intervals. Please consult `?ProfileLikelihoodProblem` on the details of the interface.

```@example example-1
using LikelihoodProfiler, Plots

# optimal values of the parameters
optpars = sol.u

# profile likelihood problem
plprob = ProfileLikelihoodProblem(optprob, optpars; profile_lower = -10., profile_upper=10., threshold = 4.0)
```

### Profile likelihood methods

LikelihoodProfiler provides a range of methods to profile likelihood functions and explore practical identifiability. The most common and simple "profiler" is the `OptimizationProfiler` method. It is based on stepwise re-optimization of the likelihood function with the constraint on the parameter (or function) of interest. We define the method and run the `solve` procedure. Please consult `?solve` on the details of the interface.

```@example example-1
method = OptimizationProfiler(optimizer = LBFGSB(), stepper = FixedStep(; initial_step=0.15))
sol = solve(plprob, method)
plot(sol, size=(800,300), margins=5Plots.mm)
```

The same `solve` interface can be used with other profiling methods. For example, a more advanced way to compute profiles is proposed by `IntegrationProfiler`. It obtains the profiles as solutions to the differential equation system. To solve this internally generated system, we need to provide a differential equations solver (`integrator`).

```@example example-1
using OrdinaryDiffEq

method = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.3,), matrix_type = :hessian)
sol = solve(plprob, method)
plot(sol, size=(800,300), margins=5Plots.mm)
```
Likelihood profiling is mostly performed to assess if the profile has intersections with the given confidence level, hence if the parameter (or function of parameters) has finite confidence interval. Another approach to the problem of practical identifiability is to compute these intersections (endpoints of the confidence interval (CI)) without restoring the full shape of the profile. One of such methods is implemented in `CICOProfiler`. It estimates CI endpoints with an optimization procedure without following the exact trajectory of the profile. 

```@example example-1
using CICOBase

method = CICOProfiler(optimizer = :LN_NELDERMEAD, scan_tol = 1e-4)
sol = solve(plprob, method)
plot(sol, size=(800,300), margins=5Plots.mm)
```