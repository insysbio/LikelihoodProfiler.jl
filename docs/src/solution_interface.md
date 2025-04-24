## Solution interface

`PLSolution` type is designed to contain the results of a profile likelihood analysis.

```@docs; canonical=false
LikelihoodProfiler.PLSolution 
```

### Retcodes

`sol::PLSolution` outputs the following retcodes, which are accessible with `get_retcodes(sol[i])` function:

- `:Identifiable` - the profile has intersection with the predefined `threshold`.
- `:NonIdentifiable` - the profile doesn't intersect the predefined `threshold`.
- `:MaxIters` - maximum number of iterations reached while computing the profile. See `maxiters` argument to the [`profile`](@ref).
- `:Failure` - the solver (optimizer or integrator) reported failure status, profiling was interrupted. 

### Visualization and tabular representation

The recipes are defined to visualize profiles saved in `sol::PLSolution` with `Plots.jl` package: `plot(sol)`, `plot(sol[i])`. 
The following keyword arguments can be used in `plot` function:

- `steps::Bool` - whether to scatter steps performed by the profiler. Defaults to `true`.
- `threshold::Bool` - whether to plot `threshold` defined in `PLProblem`. Defaults to `isfinite(threshold)`

Also each profile contained in the `sol::PLSolution` can be represented as a DataFrame with `DataFrame(sol[i])`.
