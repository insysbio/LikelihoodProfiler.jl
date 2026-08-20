
################################################## OPTIMIZATION PROFILER STEP ##################################################

function profiler_step!(profiler_cache::ProfilerCache, target::ParameterTarget, solver_cache::OptimizationSolverCache)
  @unpack opt_cache, stepper = solver_cache
  idx = get_profile_idx(profiler_cache)

  pars_guess = propose_next_pars!(profiler_cache, stepper)

  sol = solve_profile_proposal!(opt_cache, pars_guess, idx)

  if !SciMLBase.successful_retcode(sol.retcode) && stepper isa AdaptiveStep
    pars_cur = profiler_cache.pars_cur
    failed_x = pars_guess[idx]
    @. profiler_cache.pars_cache = pars_cur + (pars_guess - pars_cur) / 2
    pars_guess = clamp_profile_proposal!(profiler_cache.pars_cache, profiler_cache)

    if pars_guess[idx] != pars_cur[idx]
      @warn "Optimization failed at profile point x = $failed_x. Retrying once with half the profile step."
      sol = solve_profile_proposal!(opt_cache, pars_guess, idx)
    end
  end

  if SciMLBase.successful_retcode(sol.retcode)
    fill_x_full!(profiler_cache.pars_cur, sol.u, idx, pars_guess[idx])
    profiler_cache.x_cur = pars_guess[idx]
    profiler_cache.obj_cur = sol.objective
    profiler_cache.iter += 1
  else
    @warn "Solver returned $(sol.retcode) retcode at proposed profile point x = $(pars_guess[idx]). Profiling is interrupted."
  end
  solver_cache.retcode_original = sol.retcode
  return nothing
end

function solve_profile_proposal!(opt_cache::OptimizationCache, pars_guess, idx)
  fill_x_reduced!(opt_cache.reinit_cache.u0, pars_guess, idx)
  set_x_fixed!(opt_cache.reinit_cache.p, pars_guess[idx])
  return solve_opt_cache(opt_cache)
end

function solve_opt_cache(opt_cache::OptimizationCache)
  if isempty(opt_cache.reinit_cache.u0)
    return solve_empty_opt_cache(opt_cache)
  else
    return solve!(opt_cache)
  end
end

function solve_empty_opt_cache(opt_cache::OptimizationCache)
  u = opt_cache.reinit_cache.u0
  p = opt_cache.reinit_cache.p

  t = @elapsed obj = opt_cache.f(u, p)
  stats = SciMLBase.OptimizationStats(; iterations = 1, time = t, fevals = 1)
  retcode = isfinite_objective(obj) ? ReturnCode.Success : ReturnCode.Failure

  return SciMLBase.build_solution(opt_cache, opt_cache.opt, u, obj;
    retcode, stats)
end

################################################## INTEGRATION PROFILER STEP ##################################################

function profiler_step!(profiler_cache::ProfilerCache, target::AbstractProfileTarget, solver_cache::IntegrationSolverCache)
  @unpack ode_cache, opt_cache = solver_cache

  integrator = ode_cache
  SciMLBase.step!(integrator)

  if SciMLBase.successful_retcode(integrator.sol.retcode)
    profiler_cache.pars_cur .= view(integrator.u, 1:length(integrator.u)-1)
    profiler_cache.x_cur = integrator.t
    profiler_cache.obj_cur = evaluate_obj(get_plprob(profiler_cache), profiler_cache.pars_cur)
    profiler_cache.iter += 1
  else
    @warn "Solver returned $(integrator.sol.retcode) retcode at profile point x = $(profiler_cache.x_cur). Profiling is interrupted."

  end
  solver_cache.retcode_original = integrator.sol.retcode
  return nothing
end

################################################## OPTIMIZATION PROFILER STEPPERS ##################################################

"""
    AdaptiveInitialStep(; rel_step=0.01, abs_step=1e-3, min_step=1e-4, max_step=Inf)

Initial profile step rule that adapts to the current profiled parameter value.

The proposed step is `max(rel_step * abs(x), abs_step)`, where `x` is the
current profiled parameter value. The result is finally clamped to
`[min_step, max_step]`.

This is useful when parameters have different numerical scales and a single fixed
absolute step would be inefficient.
"""
struct AdaptiveInitialStep{T}
    rel_step::T
    abs_step::T
    min_step::T
    max_step::T
end

function AdaptiveInitialStep(; rel_step=0.01, abs_step=1e-3,
                             min_step=1e-4, max_step=Inf)
    rel_step > 0 || throw(ArgumentError("rel_step must be positive."))
    abs_step > 0 || throw(ArgumentError("abs_step must be positive."))
    min_step > 0 || throw(ArgumentError("min_step must be positive."))
    max_step >= min_step || throw(ArgumentError("max_step must be >= min_step."))

    T = float(promote_type(typeof(rel_step), typeof(abs_step),
                           typeof(min_step), typeof(max_step)))
    return AdaptiveInitialStep(T(rel_step), T(abs_step), T(min_step), T(max_step))
end


"""
    FixedStep{S}

Profiler stepper that proposes each profile point from the configured `initial_step`
rule without adapting to trial objective values.

### Constructors

```julia
FixedStep(;initial_step=AdaptiveInitialStep())
```

### Keyword arguments
- `initial_step=AdaptiveInitialStep()`: The step rule to use for each profile step. This can be a number (for a constant absolute step size), an `AdaptiveInitialStep`, or a callable `ctx -> step` for custom logic depending on the current profiler cache.
  If a number is provided, it is automatically wrapped as a function.
"""
struct FixedStep{S} <: AbstractProfilerStep
  initial_step::S
end

FixedStep(; initial_step=AdaptiveInitialStep()) = FixedStep(prepare_initial_step(initial_step))

"""
    AdaptiveStep(; initial_step=AdaptiveInitialStep(),
                   predictor=LinearPredictor(),
                   controller=ObjectiveStepControl())

Profiler stepper that adapts the profile step length based on the observed objective
increase from trial points.

The first profile step uses `initial_step`. Later steps use `predictor` to choose the
direction in parameter space and `controller` to keep the next objective increase in
a useful range.
"""
struct AdaptiveStep{S,P,C} <: AbstractProfilerStep
    initial_step::S
    predictor::P
    controller::C
end

"""
    ObjectiveStepControl(; threshold_fraction=0.1,
                           target_factor=1.5,
                           lower_factor=0.25,
                           min_obj_step=1e-2,
                           max_obj_step=Inf,
                           min_x_step=1e-4,
                           max_x_step=Inf,
                           step_factor=1.5,
                           max_x_step_growth=5.0,
                           maxiters=15)

Controls adaptive profile stepping by defining acceptable objective and profile-axis
step ranges.

For finite likelihood thresholds, the target objective increase is based on
`threshold_fraction * threshold`. For infinite thresholds, the target
increase is based on `target_factor * abs(obj_cur - obj_prev)`.

The previous optimized profile objective change is used to nudge the first
trial step up or down. Trial points below the upper objective target are then
treated as safe to grow until the upper target, profile bound, or iteration
limit is reached. The next profile-axis step is additionally capped by
`max_x_step_growth * previous_step` to avoid abrupt jumps in flat profile regions.
"""
struct ObjectiveStepControl{T}
  threshold_fraction::T
  target_factor::T
  lower_factor::T

  min_obj_step::T
  max_obj_step::T

  min_x_step::T
  max_x_step::T

  step_factor::T
  max_x_step_growth::T
  maxiters::Int
end

function ObjectiveStepControl(;
  threshold_fraction = 0.1,
    target_factor = 1.5,
    lower_factor = 0.25,

    min_obj_step = 1e-2,
    max_obj_step = Inf,

    min_x_step = 1e-4,
    max_x_step = Inf,
    max_x_step_growth = 5.0,

    step_factor = 1.5,
    maxiters = 15,
)

  0 < threshold_fraction <= 1 || throw(ArgumentError("Threshold fraction must be between 0 and 1."))
  target_factor > 0 || throw(ArgumentError("Target factor must be positive."))
  0 <= lower_factor <= 1 || throw(ArgumentError("Lower factor must be between 0 and 1."))
  min_obj_step >= 0 || throw(ArgumentError("Minimum objective step must be non-negative."))
  max_obj_step >= min_obj_step || throw(ArgumentError("Maximum objective step must be greater than or equal to minimum objective step."))
  min_x_step > 0 || throw(ArgumentError("Minimum x step must be positive."))
  max_x_step >= min_x_step || throw(ArgumentError("Maximum x step must be greater than or equal to minimum x step."))
  step_factor > 1 || throw(ArgumentError("Step factor must be greater than 1."))
  max_x_step_growth >= 1 || throw(ArgumentError("Maximum step growth must be greater than or equal to 1."))
  maxiters > 0 || throw(ArgumentError("Maximum iterations must be positive."))

  T = float(promote_type(
    typeof(threshold_fraction), typeof(target_factor), typeof(lower_factor),
    typeof(min_obj_step), typeof(max_obj_step),
    typeof(min_x_step), typeof(max_x_step),
    typeof(step_factor), typeof(max_x_step_growth),
  ))

  return ObjectiveStepControl(
    T(threshold_fraction),
    T(target_factor),
    T(lower_factor),
    T(min_obj_step),
    T(max_obj_step),
    T(min_x_step),
    T(max_x_step),
    T(step_factor),
    T(max_x_step_growth),
    Int(maxiters),
  )
end

"""
    LinearPredictor()

Extrapolates all parameters using linear extrapolation based on the last two successful profile points.
"""
struct LinearPredictor <: AbstractStepPredictor end
"""
    SingleAxisPredictor()
Extrapolates only the profiled parameter using linear extrapolation based on the last two successful profile points.
"""
struct SingleAxisPredictor <: AbstractStepPredictor end

AdaptiveStep(; initial_step=AdaptiveInitialStep(), predictor=LinearPredictor(), controller=ObjectiveStepControl()) =
  AdaptiveStep(prepare_initial_step(initial_step), predictor, controller)

function prepare_initial_step(step::Number)
  step isa Real || throw(ArgumentError("Initial step size must be real."))
  isfinite(step) && step > 0 || throw(ArgumentError("Initial step size must be finite and positive."))
  value = float(step)
  return _ -> value
end
prepare_initial_step(step::Function) = step
prepare_initial_step(step::AdaptiveInitialStep) = step
prepare_initial_step(step) = ctx -> begin
  applicable(step, ctx) || throw(ArgumentError("Initial step rule must be callable with the profiler cache."))
  step(ctx)
end

function get_initial_step(step::AbstractProfilerStep, ctx)
  value = get_initial_step(step.initial_step, ctx)
  value isa Real || throw(ArgumentError("Initial step rule must return a real number."))
  isfinite(value) && value > 0 || throw(ArgumentError("Initial step rule must return a finite, positive step size."))
  return float(value)
end

function get_initial_step(s::AdaptiveInitialStep, profiler_cache::ProfilerCache)
  x = profiler_cache.x_cur
  step = max(s.rel_step * abs(x), s.abs_step)
  return clamp(step, s.min_step, s.max_step)
end

get_initial_step(f::Function, profiler_cache::ProfilerCache) = f(profiler_cache)

function propose_next_pars!(profiler_cache::ProfilerCache, s::FixedStep)
  @unpack pars_cur, x_cur, profile_range = profiler_cache
  idx = get_profile_idx(profiler_cache)
  step_size = get_initial_step(s, profiler_cache)
  dir = get_profile_dir(profiler_cache)

  profiler_cache.pars_cache .= pars_cur
  profiler_cache.pars_cache[idx] = clamp_within_bounds(x_cur + dir * step_size, profile_range)
  return profiler_cache.pars_cache
end

function propose_next_pars!(profiler_cache::ProfilerCache, s::AdaptiveStep)
  @unpack pars_cur = profiler_cache

  isfirststep = profiler_cache.iter < 2

  if isfirststep
    # If this is the first step, we can use the previous pars to compute the direction
    return propose_next_pars!(profiler_cache, FixedStep(; initial_step=s.initial_step))
  else
    step_dir = get_step_dir(profiler_cache, s.predictor)
    step_size, retcode = adapt_step_size(profiler_cache, s.controller, step_dir)
    if retcode == :Success
      return apply_step!(profiler_cache, pars_cur, step_dir, step_size)
    elseif retcode in (:Failure, :MaxIters, :MinStep) && !(s.predictor isa SingleAxisPredictor)
      @warn "Adaptive stepping returned $retcode with $(typeof(s.predictor)). Retrying with SingleAxisPredictor."
      step_dir = get_step_dir(profiler_cache, SingleAxisPredictor())
      step_size, retcode = adapt_step_size(profiler_cache, s.controller, step_dir)
      retcode == :Success && return apply_step!(profiler_cache, pars_cur, step_dir, step_size)
    end

    if retcode == :MinStep
      @warn "Adaptive stepping reached the minimum profile step. Continuing with that step."
      return apply_step!(profiler_cache, pars_cur, step_dir, step_size)
    end

    @warn "Adaptive stepping did not find a suitable step size and returned $retcode. Continuing with the previous successful step size."
    x_cur = profiler_cache.x_cur
    x_prev = get_x_prev(profiler_cache)
    step_prev = abs(x_cur - x_prev)
    return propose_next_pars!(profiler_cache, FixedStep(; initial_step=step_prev))
  end
end

function apply_step!(profiler_cache::ProfilerCache, pars_cur, step_dir, step_size)
  @. profiler_cache.pars_cache = pars_cur + step_dir * step_size
  clamp_profile_proposal!(profiler_cache.pars_cache, profiler_cache)
  return profiler_cache.pars_cache
end

function get_step_dir(profiler_cache::ProfilerCache, ::SingleAxisPredictor)
  idx = get_profile_idx(profiler_cache)
  x_dir = get_profile_dir(profiler_cache)

  pars_cur = profiler_cache.pars_cur
  step_dir = zeros(eltype(pars_cur), length(pars_cur))
  step_dir[idx] = x_dir
  return step_dir
end

function get_step_dir(profiler_cache::ProfilerCache, ::LinearPredictor)
  @unpack pars_cur, x_cur = profiler_cache
  pars_prev = get_pars_prev(profiler_cache)
  x_prev = get_x_prev(profiler_cache)
  dx = abs(x_cur - x_prev)

  if dx <= eps(max(abs(x_cur), abs(x_prev)))
    @warn "The last two profile points are too close to each other (dx = $dx). Using SingleAxisPredictor instead of LinearPredictor."
    return get_step_dir(profiler_cache, SingleAxisPredictor())
  end
  return (pars_cur .- pars_prev) ./ dx
end

function adapt_step_size(profiler_cache::ProfilerCache, step_control::ObjectiveStepControl, step_dir)
  obj_low, obj_high = objective_step_bounds(step_control, profiler_cache)
  step_prev = initial_adaptive_step_size(profiler_cache, step_control, obj_low, obj_high)

  obj = trial_step_obj!(profiler_cache, step_dir, step_prev)

  if !isfinite_objective(obj)
    @warn "Objective function returned $obj at profile point x = $(profiler_cache.x_cur)."
    return step_prev, :Failure
  end

  step_factor = obj > obj_high ? 1 / step_control.step_factor : step_control.step_factor
  obj_prev = obj

  for _ in 1:step_control.maxiters
    step_next = clamp_x_step(profiler_cache, step_control, step_prev * step_factor)
    if step_next == step_prev
      retcode = step_factor > 1 ? :Success : :MinStep
      return step_prev, retcode
    end

    obj_next = trial_step_obj!(profiler_cache, step_dir, step_next)

    if !isfinite_objective(obj_next)
      @warn "Objective function returned $obj_next at profile point x = $(profiler_cache.x_cur)."
      return step_next, :Failure
    elseif step_factor > 1 && obj_next >= obj_high
      return interpolate_step_size(step_prev, step_next, obj_prev, obj_next, obj_high), :Success
    elseif step_factor < 1 && obj_next <= obj_high
      return interpolate_step_size(step_prev, step_next, obj_prev, obj_next, obj_high), :Success
    end

    step_prev = step_next
    obj_prev = obj_next
  end

  if step_factor > 1
    return step_prev, :Success # largest evaluated step
  else
    return step_prev, :MaxIters
  end
end

function initial_adaptive_step_size(profiler_cache::ProfilerCache, step_control::ObjectiveStepControl, obj_low, obj_high)
  step = get_step_prev(profiler_cache)
  obj_cur = profiler_cache.obj_cur
  obj_step_prev = obj_cur - get_obj_prev(profiler_cache)
  lower = obj_low - obj_cur
  upper = obj_high - obj_cur

  if obj_step_prev > upper
    step /= step_control.step_factor
  elseif -lower <= obj_step_prev < lower
    step *= step_control.step_factor
  end

  return clamp_x_step(profiler_cache, step_control, step)
end

function clamp_x_step(profiler_cache::ProfilerCache, step_control::ObjectiveStepControl, step)
  max_x_step = get_max_x_step(profiler_cache, step_control)
  min_x_step = min(step_control.min_x_step, max_x_step)
  return clamp(step, min_x_step, max_x_step)
end

function trial_step_obj!(profiler_cache::ProfilerCache, step_dir, step)
  plprob = get_plprob(profiler_cache)
  pars_cur = profiler_cache.pars_cur
  @. profiler_cache.pars_cache = pars_cur + step_dir * step
  pars = clamp_profile_proposal!(profiler_cache.pars_cache, profiler_cache)
  return evaluate_obj(plprob, pars)
end

isfinite_objective(obj) = obj isa Real && isfinite(obj)

function objective_step_bounds(step_control::ObjectiveStepControl, profiler_cache::ProfilerCache)
    threshold = get_threshold(profiler_cache)
    obj_cur = profiler_cache.obj_cur

    if isfinite(threshold)
      upper = clamp(step_control.threshold_fraction * threshold, step_control.min_obj_step, step_control.max_obj_step)
    else
      obj_prev = get_obj_prev(profiler_cache)
      upper = step_control.target_factor * abs(obj_cur - obj_prev)
      upper = clamp(upper, step_control.min_obj_step, step_control.max_obj_step)
    end

    lower = step_control.lower_factor * upper

    return obj_cur + lower, obj_cur + upper
end

function get_max_x_step(profiler_cache::ProfilerCache, step_control::ObjectiveStepControl)
  profile_lb, profile_ub = profiler_cache.profile_range
  x_cur = profiler_cache.x_cur
  bound_step = !isleft(profiler_cache) ? profile_ub - x_cur : x_cur - profile_lb
  previous_step = get_step_prev(profiler_cache)
  growth_limited_step = step_control.max_x_step_growth * previous_step
  max_x_step = min(bound_step, step_control.max_x_step, growth_limited_step)
  max(max_x_step, zero(max_x_step))
end

function get_step_prev(profiler_cache)
  x_cur = profiler_cache.x_cur
  x_prev = get_x_prev(profiler_cache)
  return abs(x_cur - x_prev)
end

# Clamp a number within specified bounds, handling various cases for bounds
clamp_within_bounds(x::Number, bounds::Tuple) = clamp_within_bounds(x, bounds[1], bounds[2])
clamp_within_bounds(x::Number, bounds::Nothing) = x
clamp_within_bounds(x::Number, lb::Nothing, ub::Nothing) = x
clamp_within_bounds(x::Number, lb::Number, ub::Number) = clamp(x, lb, ub)
clamp_within_bounds(x::Number, lb::Nothing, ub::Number) = clamp(x, -Inf, ub)
clamp_within_bounds(x::Number, lb::Number, ub::Nothing) = clamp(x, lb, Inf)

function clamp_profile_proposal!(pars::AbstractVector, profiler_cache::ProfilerCache)
  optprob = get_optprob(profiler_cache)
  clamp_within_bounds!(pars, optprob.lb, optprob.ub)
  idx = get_profile_idx(profiler_cache)
  pars[idx] = clamp_within_bounds(pars[idx], profiler_cache.profile_range)
  return pars
end

# In-place version for θ::Vector
function clamp_within_bounds!(θ::AbstractVector, bounds::Tuple)
    for i in eachindex(θ)
        θ[i] = clamp_within_bounds(θ[i], bounds)
    end
    return θ
end

function clamp_within_bounds!(θ::AbstractVector, bounds::AbstractVector)
    for i in eachindex(θ)
        θ[i] = clamp_within_bounds(θ[i], bounds[i])
    end
    return θ
end

function clamp_within_bounds!(θ::AbstractVector, lb, ub)
    for i in eachindex(θ)
        lb_i = isnothing(lb) ? nothing : lb[i]
        ub_i = isnothing(ub) ? nothing : ub[i]
        θ[i] = clamp_within_bounds(θ[i], lb_i, ub_i)
    end
    return θ
end

# Non-inplace versions for backward compatibility
clamp_within_bounds(θ::AbstractVector, bounds::Tuple) = clamp_within_bounds!(copy(θ), bounds)
clamp_within_bounds(θ::AbstractVector, bounds::AbstractVector) = clamp_within_bounds!(copy(θ), bounds)

function interpolate_step_size(step_low, step_high, obj_low, obj_high, obj_target)
  Δobj = obj_high - obj_low
  abs(Δobj) < eps(max(abs(obj_high), abs(obj_low))) && return (step_low + step_high) / 2
  t = (obj_target - obj_low) / Δobj
  step = muladd(t, step_high - step_low, step_low)
  step_min, step_max = minmax(step_low, step_high)
  return clamp(step, step_min, step_max)
end
