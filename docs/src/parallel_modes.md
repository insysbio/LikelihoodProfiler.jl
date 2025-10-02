# Parallel execution

`LikelihoodProfiler.jl` supports parallel computation to accelerate profile likelihood calculations. All implemented methods support independent computations **across parameters and profile branches ("left" and "right")**. Depending on your system configuration and workload, you can leverage either **multithreading** or **distributed processing (workers)** for parallel execution.

## Available parallel modes

Parallelization is controlled via the `parallel_type` keyword argument in the `solve` function.

```julia
solve(plprob, method; parallel_type = :none)
```

Supported values:

- `:none` (default): Run all profiling computations sequentially.
- [`:threads`](https://docs.julialang.org/en/v1/manual/multi-threading/): Use Julia's multithreading (i.e. `Base.Threads.@threads`) across available CPU threads. To enable multithreading, set the `JULIA_NUM_THREADS` environment variable.
- [`:distributed`](https://docs.julialang.org/en/v1/stdlib/Distributed/): Use `Distributed.jl` to run tasks across multiple Julia processes (e.g., started via `addprocs()` or `julia -p N`).

## Distributed workers example

The following example illustrates how to compute profile likelihoods using distributed parallelism across Julia processes in a multi-core (or multi-node) environment.
```julia
using Distributed

addprocs(2)

@everywhere using LikelihoodProfiler, Optimization, ForwardDiff

@everywhere rosenbrock(x,p) = (1.0 - x[1])^2 + 100.0*(x[2] - x[1]^2)^2

x0 = [1., 1.]
optf = OptimizationFunction(rosenbrock, AutoForwardDiff())
optprob = OptimizationProblem(optf, x0)

plprob = ProfileLikelihoodProblem(optprob, x0; profile_lower=-5.0, profile_upper=5.0, threshold = 1.0)
meth = OptimizationProfiler(optimizer = Optimization.LBFGS(), stepper = FixedStep(; initial_step=0.1))

sol = solve(plprob, meth; parallel_type=:distributed)
```