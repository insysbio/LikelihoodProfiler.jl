
################################################## OPTIMIZATION PROFILER STEP ##################################################

function profiler_step!(profiler_cache::ProfilerCache, target::ParameterTarget, solver_cache::OptimizationSolverCache)
  @unpack opt_cache, stepper = solver_cache
  idx = get_profile_idx(profiler_cache)

  pars_guess = propose_next_pars!(profiler_cache, stepper)

  fill_x_reduced!(opt_cache.reinit_cache.u0, pars_guess, idx)
  set_x_fixed!(opt_cache.reinit_cache.p, pars_guess[idx])
  
  sol = solve_opt_cache(opt_cache)

  if SciMLBase.successful_retcode(sol.retcode)
    fill_x_full!(profiler_cache.pars_cur, sol.u, idx, pars_guess[idx])
    profiler_cache.x_cur = pars_guess[idx]
    profiler_cache.obj_cur = sol.objective
    profiler_cache.iter += 1
  else
    @warn "Solver returned $(sol.retcode) retcode at profile point x = $(profiler_cache.x_cur). Profiling is interrupted."
  end
  solver_cache.retcode_original = sol.retcode
  return nothing
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

  return SciMLBase.build_solution(opt_cache, opt_cache.opt, u, obj; 
    retcode = ReturnCode.Success, stats)
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

const DEFAULT_INIT_STEP = 1e-2

"""
    FixedStep{S}

Profiler stepper that always proposes a fixed step size in the profiling direction.

### Constructors

```julia
FixedStep(;initial_step=DEFAULT_INIT_STEP)
```

### Keyword arguments
- `initial_step=DEFAULT_INIT_STEP`: The initial (and constant) step size to use for each profile step. This can be a number (for a constant step size) or a function `(pars, idx) -> step` for custom logic depending on the current parameters and index.
  If a number is provided, it is automatically wrapped as a function.
"""
struct FixedStep{S<:Function} <: AbstractProfilerStep{S}
  initial_step::S
end

FixedStep(; initial_step=DEFAULT_INIT_STEP) = FixedStep(prepare_initial_step(initial_step))


"""
    LineSearchStep{S, L}

Profiler stepper that uses a line search to adaptively determine the step size in the direction which is chosen by the `direction` keyword argument.

### Constructors

```julia
LineSearchStep(;initial_step=DEFAULT_INIT_STEP, direction=:Secant, linesearch=InterpolationLineSearch())
```

### Keyword arguments
- `initial_step=DEFAULT_INIT_STEP`: Initial guess for the step size. Can be a number (for a constant guess) or a function `(pars, idx) -> step` for custom logic depending on the current parameters and index. If a number is provided, it is automatically wrapped as a function.
- `direction::Symbol=:Secant`: Strategy for choosing the direction of the next step. Options:
    - `:SingleAxis`: Move along the current profiling parameter only.
    - `:Secant`: Use the secant direction, i.e., the line connecting the last two points in parameter space (Default).
    - `:Gradient`: Use the gradient of the objective function as the direction.
  The choice affects how the next step is proposed in parameter space.
- `linesearch::L=InterpolationLineSearch()`: Line search algorithm (e.g., `InterpolationLineSearch()`).
"""
struct LineSearchStep{S<:Function, L} <: AbstractProfilerStep{S}
  initial_step::S
  direction::Symbol
  linesearch::L
end

function LineSearchStep(; initial_step=DEFAULT_INIT_STEP, direction=:Secant, linesearch=InterpolationLineSearch()) 
  @assert direction in (:SingleAxis, :Secant, :Gradient) "Invalid direction: $direction. Must be one of (:SingleAxis, :Secant, :Gradient)."
  LineSearchStep(prepare_initial_step(initial_step), direction, linesearch)
end


"""
    InterpolationLineSearch

Interpolation-based line search algorithm used in profile stepping.

### Constructors

```julia
InterpolationLineSearch(; objective_factor=1.25, step_size_factor=2.0, maxiters=10) 
```

### Keyword arguments
- `objective_factor::Float64 = 1.25`: Factor by which the change in the objective function from the last step is increased/decreased to set the target for the next step.
- `step_size_factor::Float64 = 2.0`: Multiplicative factor for increasing or decreasing the step size during the search.
- `maxiters::Int = 10`: Maximum number of line search iterations allowed.
"""
struct InterpolationLineSearch
  objective_factor::Float64
  step_size_factor::Float64
  maxiters::Int
end

function InterpolationLineSearch(; objective_factor=1.25, step_size_factor=2.0, maxiters=10) 
  @assert objective_factor > 1 "Objective scaling factor must be greater than 1."
  @assert step_size_factor > 0 "Step size factor must be positive."
  @assert maxiters > 0 "Maximum iterations must be positive."
  
  return InterpolationLineSearch(objective_factor, step_size_factor, maxiters)
end

#=
# TODO 
struct TrustRegionStep{S} <: AbstractProfilerStep{S}
  initial_step::S
  min_x_step::Float64
  max_x_step::Float64
  min_obj_step::Float64
  max_obj_step::Float64
end
=#
function prepare_initial_step(step::Number)
  @assert step > 0 "Initial step size must be positive."
  return (p0,i) -> float(step)
end
prepare_initial_step(step::Function) = step

get_step(s::AbstractProfilerStep{S}, pars, i) where S <: Function = s.initial_step(pars, i)


function propose_next_pars!(profiler_cache::ProfilerCache, s::FixedStep)
  @unpack pars_cur, x_cur, profile_range = profiler_cache
  idx = get_profile_idx(profiler_cache)
  step_size = get_step(s, pars_cur, idx)
  dir = get_profile_dir(profiler_cache)

  profiler_cache.pars_cache .= pars_cur
  profiler_cache.pars_cache[idx] = clamp_within_bounds(x_cur + dir * step_size, profile_range)
  return profiler_cache.pars_cache
end

function propose_next_pars!(profiler_cache::ProfilerCache, s::LineSearchStep)
  @unpack pars_cur, x_cur, profile_range = profiler_cache

  isfirststep = profiler_cache.iter < 2

  if isfirststep
    # If this is the first step, we can use the previous pars to compute the direction
    return propose_next_pars!(profiler_cache, FixedStep(; initial_step=s.initial_step))
  else
    ls_dir = compute_direction(profiler_cache, Val(s.direction))
    step_size, ls_retcode = compute_step_size(profiler_cache, s.linesearch, ls_dir)
    #@show step_size, ls_retcode
    if ls_retcode == :Success
      @. profiler_cache.pars_cache = pars_cur + ls_dir*step_size
      clamp_within_bounds!(profiler_cache.pars_cache, profile_range)
      return profiler_cache.pars_cache
    else
      @warn "Line search didn't find a suitable step size and returned $ls_retcode. Continuing with the previous successful step size (FixedStep)."
      x_cur = profiler_cache.x_cur
      x_prev = get_x_prev(profiler_cache)
      step_prev = abs(x_cur - x_prev)
      return propose_next_pars!(profiler_cache, FixedStep(; initial_step=step_prev))
    end
  end
end

function compute_direction(profiler_cache::ProfilerCache, direction::Val{:SingleAxis})
  idx = get_profile_idx(profiler_cache)
  x_dir = get_profile_dir(profiler_cache)

  pars_cur = profiler_cache.pars_cur
  ls_dir = zeros(eltype(pars_cur), length(pars_cur))
  ls_dir[idx] = x_dir
  return ls_dir
end

function compute_direction(profiler_cache::ProfilerCache, direction::Val{:Secant})
  @unpack pars_cur, x_cur = profiler_cache
  θ_prev = get_pars_prev(profiler_cache)
  x_prev = get_x_prev(profiler_cache)

  return (pars_cur .- θ_prev) ./ abs(x_cur - x_prev)
end

function compute_direction(profiler_cache::ProfilerCache, direction::Val{:Gradient})
  optprob = get_optprob(profiler_cache)
  θ_idx = get_profile_idx(profiler_cache)
  pars_cur = profiler_cache.pars_cur

  grad = evaluate_gradf(optprob, pars_cur)
  return grad ./ abs(grad[θ_idx])
end

function compute_step_size(profiler_cache::ProfilerCache, ls_alg::InterpolationLineSearch, ls_dir)
  @unpack pars_cur, x_cur, obj_cur, profile_range = profiler_cache
  x_prev = get_x_prev(profiler_cache)
  plprob = get_plprob(profiler_cache)
  obj_prev = get_obj_prev(profiler_cache)
  Δ_obj_prev = abs(obj_cur - obj_prev)

  # try previous step_size (extrapolate_profile)
  step_next = abs(x_cur - x_prev)
  @. profiler_cache.pars_cache = pars_cur + ls_dir*step_next
  θ_next = clamp_within_bounds!(profiler_cache.pars_cache, profile_range)
  obj_next = evaluate_obj(plprob, profiler_cache.pars_cache)

  increasing_profile = obj_next > obj_cur
  obj_target = increasing_profile ? obj_cur + ls_alg.objective_factor*Δ_obj_prev : obj_cur - ls_alg.objective_factor*Δ_obj_prev
  #obj_target = increasing_profile ? obj_cur + ls_alg.objective_factor*obj_cur : obj_cur - ls_alg.objective_factor*obj_cur

  # TODO add tolerance arg to the alg
  (abs(obj_next - obj_target) < 1e-6) && return step_next, :Success
  
  undershoot = increasing_profile ? 
    (obj_val -> obj_val < obj_target) : 
    (obj_val -> obj_val > obj_target)

  bracketing_cond = increasing_profile ?
    ((obj_low, obj_high) -> obj_low < obj_target < obj_high) :
    ((obj_low, obj_high) -> obj_low > obj_target > obj_high)

  #increase_step = undershoot(obj_next)
  #factor = increase_step ? ls_alg.step_size_factor : 1/ls_alg.step_size_factor
  factor = ls_alg.step_size_factor
  iter = 0
  # tmp fixing step_low, obj_low to current point
  step_low, obj_low = 0.0, obj_cur
  step_high, obj_high = step_next, obj_next
  while !bracketing_cond(obj_low, obj_high)
    iter += 1
    if iter > ls_alg.maxiters
      #@show step_low, step_high, obj_low, obj_high, obj_target
      return step_low, :MaxIters
    end
    
    step_next *= factor
    @. profiler_cache.pars_cache = pars_cur + ls_dir * step_next
    θ_next = clamp_within_bounds!(profiler_cache.pars_cache, profile_range)
    obj_next = evaluate_obj(plprob, θ_next)
    #@show step_next, obj_next, obj_target
    
    # TODO add tolerance arg to the alg
    (abs(obj_next - obj_target) < 1e-6) && return step_next, :Success

    # tmp removing undershoot condition
    step_high, obj_high = step_next, obj_next
    #=
    if undershoot(obj_next)
      step_low, obj_low = step_next, obj_next
    else
      step_high, obj_high = step_next, obj_next
    end
    =#
  end
  return interpolate_step_size(step_low, step_high, obj_low, obj_high, obj_target), :Success
end

#=
function compute_step_size(profiler_cache::ProfilerCache, ls_alg::InterpolationLineSearch, ls_dir)
  pars_cur = get_pars_cur(profiler_cache)
  x_cur = get_x_cur(profiler_cache)
  x_prev = get_x_prev(profiler_cache)
  profile_range = get_profile_range(profiler_cache)
  opt_prob = get_optprob(profiler_cache)
  obj_cur = get_obj_cur(profiler_cache)
  obj_prev = get_obj_prev(profiler_cache)
  Δ_obj_prev = abs(obj_cur - obj_prev)

  # try previous step_size (extrapolate_profile)
  step_next = abs(x_cur - x_prev)
  @. profiler_cache.pars_cache = pars_cur + ls_dir*step_next
  θ_next = clamp_within_bounds!(profiler_cache.pars_cache, profile_range)
  obj_next = evaluate_obj(opt_prob, profiler_cache.pars_cache)

  increasing_profile = obj_next > obj_cur
  obj_target = increasing_profile ? obj_cur + ls_alg.objective_factor*Δ_obj_prev : obj_cur - ls_alg.objective_factor*Δ_obj_prev
  #obj_target = increasing_profile ? obj_cur + ls_alg.objective_factor*obj_cur : obj_cur - ls_alg.objective_factor*obj_cur

  # TODO add tolerance arg to the alg
  (abs(obj_next - obj_target) < 1e-6) && return step_next, :Success
  
  undershoot = increasing_profile ? 
    (obj_val -> obj_val < obj_target) : 
    (obj_val -> obj_val > obj_target)

  bracketing_cond = increasing_profile ?
    ((obj_low, obj_high) -> obj_low < obj_target < obj_high) :
    ((obj_low, obj_high) -> obj_low > obj_target > obj_high)

  increase_step = undershoot(obj_next)
  factor = increase_step ? ls_alg.step_size_factor : 1/ls_alg.step_size_factor

  iter = 0
  step_low, obj_low = step_next, obj_next
  step_high, obj_high = step_next, obj_next
  while !bracketing_cond(obj_low, obj_high)
    iter += 1
    if iter > ls_alg.maxiters
      @show step_low, step_high, obj_low, obj_high, obj_target
      return step_low, :MaxIters
    end
    
    step_next *= factor
    @. profiler_cache.pars_cache = pars_cur + ls_dir * step_next
    θ_next = clamp_within_bounds!(profiler_cache.pars_cache, profile_range)
    obj_next = evaluate_obj(opt_prob, θ_next)
    @show step_next, obj_next, obj_target
    # TODO add tolerance arg to the alg
    (abs(obj_next - obj_target) < 1e-6) && return step_next, :Success

    if undershoot(obj_next)
      step_low, obj_low = step_next, obj_next
    else
      step_high, obj_high = step_next, obj_next
    end
    
  end
  return interpolate_step_size(step_low, step_high, obj_low, obj_high, obj_target), :Success
end
=#

# Clamp a number within specified bounds, handling various cases for bounds
clamp_within_bounds(x::Number, bounds::Tuple) = clamp_within_bounds(x, bounds[1], bounds[2])
clamp_within_bounds(x::Number, bounds::Nothing) = x
clamp_within_bounds(x::Number, lb::Nothing, ub::Nothing) = x
clamp_within_bounds(x::Number, lb::Number, ub::Number) = clamp(x, lb, ub)
clamp_within_bounds(x::Number, lb::Nothing, ub::Number) = clamp(x, -Inf, ub)
clamp_within_bounds(x::Number, lb::Number, ub::Nothing) = clamp(x, lb, Inf)

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

# Non-inplace versions for backward compatibility
clamp_within_bounds(θ::AbstractVector, bounds::Tuple) = clamp_within_bounds!(copy(θ), bounds)
clamp_within_bounds(θ::AbstractVector, bounds::AbstractVector) = clamp_within_bounds!(copy(θ), bounds)

function interpolate_step_size(step_low, step_high, obj_low, obj_high, obj_target)
  if abs(obj_high - obj_low) < eps()
    return (step_low + step_high) / 2
  end
  step_interp = step_low + (obj_target - obj_low) / (obj_high - obj_low) * (step_high - step_low)
  return clamp(step_interp, step_low, step_high)
end