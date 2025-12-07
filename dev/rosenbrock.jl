using Plots

rosenbrock(x,p) = (1.0 - x[1])^2 + 100.0*(x[2] - x[1]^2)^2

x = range(-2.0, 2.0, length=400)
y = range(-1.0, 3.0, length=400)
z = [rosenbrock([x,y], nothing) for y in y, x in x]
p = contourf(x, y, z; levels=30, c = :thermal,  colorbar = true, title="Rosenbrock Contours",)

f1(x,p) = x[2] + x[1]^2
f2(x,p) = x[2] - x[1]^2             # Rosenbrock valley
f3(x,p) = x[2] - x[1]^2 + 0.2*x[1]  # tilted valley

# Line for the Rosenbrock valley: x2 = x1^2
x1_line = range(-2, 2, length=400)
valley = x1_line .^ 2

# Curves for f1, f2, f3 for visualization
# (We plot level sets passing through the optimum at x* = (1, 1))
x1_curve = x1_line
f1_curve = f1.([(x,y) for (x,y) in zip(x,y)], nothing)
f2_curve = f2.([x for x in zip(x1_curve, x1_curve.^2)], nothing)
f3_curve = f3.([x for x in zip(x1_curve, x1_curve.^2)], nothing)

plot!(p, x1_curve, f1_curve, label="f1", ylims = (-1, 3), color=:blue, lw=2, linestyle=:dash)
plot!(p, x1_curve, f2_curve, label="f2 (valley)", color=:green, lw=2, linestyle=:dash)


#=

### Profiling Functions of Parameters

In addition to profiling individual parameters, `LikelihoodProfiler` can profile arbitrary scalar functions of the parameter vector. This is useful for studying:
- identifiable and non-identifiable combinations of parameters,
- reparameterizations,
- predictions (predictability analysis) or derived-quantity identifiability,
- curvature and valley structure of the objective function.

To illustrate this functionality, we consider three simple functions of the Rosenbrock parameters. 
1. `f1(x,p) = x[2] + x[1]^2` 
2. `f2(x,p) = x[2] - x[1]^2` 
3. `f3(x,p) = x[2] - x[1]^2 + 0.2*x[1]` 

```julia
f1(x,p) = x[2] + x[1]^2
f2(x,p) = x[2] - x[1]^2
f3(x,p) = x[2] - x[1]^2 + 0.2*x[1]

optf1 = OptimizationFunction(f1, AutoForwardDiff())
optf2 = OptimizationFunction(f2, AutoForwardDiff())
optf3 = OptimizationFunction(f3, AutoForwardDiff())

func_plprob = ProfileLikelihoodProblem(optprob, optpars, [optf1, optf2, optf3];
  profile_lower=-10.0, profile_upper=10.0)

func_meth = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.03,), matrix_type = :hessian)

sol_func = solve(func_plprob, func_meth)
```
=#