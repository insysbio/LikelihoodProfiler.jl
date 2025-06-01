# Parallel execution

`LikelihoodProfiler.jl` supports parallel computation to accelerate profile likelihood calculations. All implemented methods support independent computations across parameters and profile branches ("left" and "right"). Depending on your system configuration and workload, you can leverage either **multithreading** or **distributed processing (workers)** for parallel execution.

## Available parallel modes

Parallelization is controlled via the `parallel_type` keyword argument in the `profile` function.

```julia
profile(plprob, method; parallel_type = :none)
```

Supported values:

- `:none` (default): Run all profiling computations sequentially.
- `:threads`: Use Julia's multithreading (i.e. `Base.Threads.@threads`) across available CPU threads.
- `:distributed`: Use `Distributed.jl` to run tasks across multiple Julia processes (e.g., started via `addprocs()` or `julia -p N`).