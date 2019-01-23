# Visualization

[`LikelihoodProfiler.get_interval`](@ref) function outputs estimated
confidence interval along with other data as
[`LikelihoodProfiler.ParamInterval`](@ref) structure.

**LikelihoodProfiler** provides a `@recipe` for **Plots.jl** to visualize
confidence interval estimation and plot parameter profile based on
[`LikelihoodProfiler.ParamInterval`](@ref).

```
using LikelihoodProfiler

# Likelihood function
f(x) = 5.0 + (x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2

# Calculate parameters intervals for x[1], x[2], x[3]
res = [
    get_interval(
        [3., 2., 2.1],
        i,
        f,
        :CICO_ONE_PASS;
        loss_crit = 9.
    ) for i in 1:3]

# Plot parameter profile x[1]
using Plots
plotly()
plot(res[1])
```
