
## [Profile Likelihood Methods](@id profile_likelihood_methods)

LikelihoodProfiler provides a range of methods to profile likelihood functions and explore practical identifiability. The method should be provided as the second argument to the [`solve`](@ref) function.

### [Optimization-based profiles](@id optimization_based_profiles)

The method computes profiles for each parameter by iteratively changing the value of the parameter and re-optimizing the likelihood function with respect to all other parameters. 

```@docs; canonical=false
OptimizationProfiler
```

#### Optimization steppers

`OptimizationProfiler` uses a `stepper` to choose the next value along the profiled parameter axis before re-optimizing all remaining parameters.

Use [`FixedStep`](@ref) when you want a predictable step rule that does not inspect
the objective function at trial points:

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

With its default settings, `AdaptiveInitialStep()` proposes
`max(0.01 * abs(x), 1e-3)` and clamps the result to `[1e-4, Inf]`. These
defaults are conservative for parameters on their original scale. For parameters
on a log or otherwise standardized scale, a fixed initial step such as `0.01` can
also be a convenient choice.

#### How adaptive stepping works

After the first profile point, `AdaptiveStep` separates the proposal into a
direction and a length:

1. The predictor estimates how all parameters should move when the profiled
   parameter changes.
2. The controller evaluates the objective at trial points along that direction,
   without re-optimizing the nuisance parameters.
3. The trial step is increased or decreased until its objective change is close
   to the upper target, or a profile bound or configured step limit is reached.
4. The nuisance parameters are re-optimized only after the proposal has been
   selected.

The trial objective therefore measures the quality of the optimization starting
point, not the final profile objective. The final objective will usually be lower
after re-optimization, and it is allowed to cross the profile threshold. A
threshold crossing terminates that profile branch in the usual way.

The default [`LinearPredictor`](@ref) extrapolates the complete parameter vector
from the two previous optimized profile points. It is usually the most efficient
choice for smooth profiles because it follows the estimated profile path through
the nuisance-parameter space. [`SingleAxisPredictor`](@ref) changes only the
profiled parameter and leaves all nuisance parameters at their current optimum.
It is more conservative and can be useful for profiles with abrupt turns, local
optima, or unreliable extrapolation.

```julia
# Default: first-order extrapolation of the complete profile path
adaptive = AdaptiveStep(; predictor = LinearPredictor())

# Conservative alternative: move only the profiled parameter
conservative = AdaptiveStep(; predictor = SingleAxisPredictor())
```

If adaptive search with `LinearPredictor` fails or reaches its minimum step, the
profiler automatically retries with `SingleAxisPredictor`. If optimization at an
accepted adaptive proposal fails, it is retried once with half the proposed step.

#### Objective step control

[`ObjectiveStepControl`](@ref) defines the desired objective change and the
allowed profile-axis step sizes. For a finite profile threshold, its upper target
is

```math
\Delta f_{\mathrm{target}} =
\operatorname{clamp}(\mathtt{threshold\_fraction}\,\tau,
                     \mathtt{min\_obj\_step},
                     \mathtt{max\_obj\_step}),
```

where `tau` is the objective threshold relative to the optimum. The default
`threshold_fraction=0.1` aims for approximately ten profile steps over that
objective range when the local profile shape is regular.

When `threshold=Inf`, no confidence threshold is available as an objective scale.
The target is instead based on `target_factor` times the absolute objective change
between the two previous optimized profile points, clamped by `min_obj_step` and
`max_obj_step`.

The most useful tuning options are:

| Option | Effect of increasing it |
|:--|:--|
| `threshold_fraction` | Fewer, larger steps for finite thresholds |
| `min_obj_step` | Larger steps when the objective is flat or `threshold=Inf` |
| `step_factor` | Faster but coarser trial-step growth and reduction |
| `max_x_step_growth` | Faster growth across flat profile regions |
| `min_x_step` | Prevents excessively dense sampling |
| `max_x_step` | Limits jumps along the profiled parameter axis |

For a profile that remains too dense in flat regions, first increase the initial
step or `min_obj_step`. If growth is still too slow, increase
`max_x_step_growth`. For unstable optimization or a sharply curved profile,
reduce `threshold_fraction` or `max_x_step`, or select
`SingleAxisPredictor()`.

```julia
controller = ObjectiveStepControl(
    threshold_fraction = 0.1,
    min_obj_step = 0.01,
    min_x_step = 1e-4,
    max_x_step = 0.5,
)

method = OptimizationProfiler(
    optimizer = LBFGSB(),
    stepper = AdaptiveStep(; controller),
)
```

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
