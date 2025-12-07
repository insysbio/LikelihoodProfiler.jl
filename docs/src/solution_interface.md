## [Solution Interface](@id solution_interface)

`ProfileLikelihoodSolution` type is designed to contain the results of a profile likelihood analysis.

```@docs; canonical=false
LikelihoodProfiler.ProfileLikelihoodSolution 
```

### Retcodes

`sol::ProfileLikelihoodSolution` outputs the following retcodes, which are accessible with `retcodes(sol[i])` function:

- `:Identifiable` - the profile has intersection with the predefined `threshold`.
- `:NonIdentifiable` - the profile doesn't intersect the predefined `threshold`.
- `:MaxIters` - maximum number of iterations reached while computing the profile. See `maxiters` argument to the [`solve`](@ref).
- `:Failure` - the solver (optimizer or integrator) reported failure status, profiling was interrupted. 

### Endpoints (confidence-interval estimates)

`endpoints(sol[i])` returns the estimated crossing points of the profile with the chosen likelihood `threshold` for the i-th profile. 

### Visualization and tabular representation

The recipes are defined to visualize profiles saved in `sol::ProfileLikelihoodSolution` with `Plots.jl` package: `plot(sol)`, `plot(sol[i])`. 
The following keyword arguments can be used in `plot` function:

- `steps::Bool` - whether to scatter steps performed by the profiler. Defaults to `true`.
- `threshold::Bool` - whether to plot `threshold` defined in `ProfileLikelihoodProblem`. Defaults to `isfinite(threshold)`

Also each profile contained in the `sol::ProfileLikelihoodSolution` can be represented as a DataFrame with `DataFrame(sol[i])`.
