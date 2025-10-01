## Defining a ProfileLikelihoodProblem

`ProfileLikelihoodProblem` type is designed to contain the necessary information to define a profile likelihood problem.

```@docs; canonical=false
ProfileLikelihoodProblem
```

### Targets (what you profile)

A **profile target** specifies the quantity whose profile you want: either one or more **parameters** by index, or one or more **functions of the parameters**. Targets carry their own profile bounds (`profile_lower`, `profile_upper`) which delimit where each profile is traced.

#### ParameterTarget

`ParameterTarget` profiles model **parameters by index**. 
While you can build `ParameterTarget` directly, most users prefer the convenience constructor on `ProfileLikelihoodProblem`, which accepts idxs and bounds and constructs the target for you. See `ProfileLikelihoodProblem` interface.

```@docs; canonical=false
ParameterTarget
```

##### Explicit construction example

```@example ex-1
using LikelihoodProfiler, Optimization

f = OptimizationFunction((θ,p)->sum(abs2, θ))
optprob = OptimizationProblem(f, [1.0, 2.0, 3.0])

pt = ParameterTarget(; idxs=[1,3],
                     profile_lower=[-5.0, -1.0],
                     profile_upper=[ 2.0,  4.0])

plprob = ProfileLikelihoodProblem(optprob, [0.0, 0.0, 0.0], pt;
                                  conf_level=0.95)  # threshold derived from χ²
```

#### FunctionTarget

`FunctionTarget` profiles **functions of the parameters**.
While you can build `FunctionTarget` directly, most users prefer the convenience constructor on `ProfileLikelihoodProblem`, which accepts functions and bounds and constructs the target for you. See `ProfileLikelihoodProblem` interface.

```@docs; canonical=false
FunctionTarget
```

##### Explicit construction example

```@example ex-2
using LikelihoodProfiler, Optimization

f = OptimizationFunction((θ,p)->sum(abs2, θ))
optprob = OptimizationProblem(f, [1.0, 2.0, 3.0])

g1 = OptimizationFunction((θ,p)->θ[1] + θ[2])
g2 = OptimizationFunction((θ,p)->θ[2] - θ[3])

ft = FunctionTarget(; fs=[g1,g2],
                    profile_lower=[-2.0, -1.0],
                    profile_upper=[ 2.0,  1.0])

plprob = ProfileLikelihoodProblem(optprob, [0.0, 0.0, 0.0], ft)
```

### Problem “sugar” constructors

For ergonomics, `ProfileLikelihoodProblem` provides keyword constructors that build the target for you.

```julia
ProfileLikelihoodProblem(::OptimizationProblem, ::AbstractVector{<:Real}; idxs, profile_lower, profile_upper)
ProfileLikelihoodProblem(::OptimizationProblem, ::AbstractVector{<:Real}, ::Union{OptimizationFunction,AbstractVector{<:OptimizationFunction}}; profile_lower, profile_upper)
```