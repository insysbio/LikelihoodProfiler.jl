
## [Profile Likelihood Methods](@id profile_likelihood_methods)

LikelihoodProfiler provides a range of methods to profile likelihood functions and explore practical identifiability. The method should be provided as the second argument to the [`solve`](@ref) function.

### [Optimization-based profiles](@id optimization_based_profiles)

The method computes profiles for each parameter by iteratively changing the value of the parameter and re-optimizing the likelihood function with respect to all other parameters. 

```@docs; canonical=false
OptimizationProfiler
```

#### Optimization steppers

`OptimizationProfiler` uses a `stepper` to choose the next value along the profiled parameter axis before re-optimizing all remaining parameters.

Use [`FixedStep`](@ref) when you want a predictable step size throughout the profile:

```julia
method = OptimizationProfiler(
	optimizer = LBFGSB(),
	stepper = FixedStep(; initial_step = 0.1),
)
```

Use [`AdaptiveStep`](@ref) when you want the profiler to adjust the step length based on the objective increase observed at trial points:

```julia
method = OptimizationProfiler(
	optimizer = LBFGSB(),
	stepper = AdaptiveStep(; initial_step = AdaptiveInitialStep()),
)
```

The default [`AdaptiveInitialStep`](@ref) scales the first step with the current profiled value and clamps it to configured minimum and maximum bounds. This is useful when parameters have different numerical scales.

```@docs; canonical=false
AdaptiveInitialStep
FixedStep
AdaptiveStep
ObjectiveStepControl
LinearPredictor
SingleAxisPredictor
```

### [Integration-based profiles](@id integration_based_profiles)

The method computes profiles for each parameter (or function of parameters) by integrating the differential equations system. 

```@docs; canonical=false
IntegrationProfiler
```

References:
1. Chen, J.-S. & Jennrich, R. I. Simple Accurate Approximation of Likelihood Profiles. Journal of Computational and Graphical Statistics 11, 714–732 (2002).
2. Chen, J.-S. & Jennrich, R. I. The Signed Root Deviance Profile and Confidence Intervals in Maximum Likelihood Analysis. Journal of the American Statistical Association 91, 993–998 (1996).

### [Confidence Intervals by Constrained Optimization (CICO)](@id cico_profiles)

The method computes intersections (endpoints of the confidence interval (CI)) of the profile with the predefined confidence level (`threshold`) without restoring the exact trajectory of the profile. Requires using [CICOBase](https://github.com/insysbio/CICOBase.jl) package.

```@docs; canonical=false
CICOProfiler 
```

References:
1. Borisov, I. & Metelkin, E. Confidence intervals by constrained optimization—An algorithm and software package for practical identifiability analysis in systems biology. PLoS Comput Biol 16, e1008495 (2020).
2. Venzon, D. J. & Moolgavkar, S. H. A Method for Computing Profile-Likelihood-Based Confidence Intervals. Applied Statistics 37, 87 (1988).


### [Quadratic approximation (FIM curvature at optimum)](@id fim_profiles)

```@docs; canonical=false
QuadraticApproxProfiler
```
