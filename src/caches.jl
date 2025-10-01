mutable struct OptimizationSolverCache{O,StepType,StatsType,R} <: AbstractSolverCache
    opt_cache::O
    stepper::StepType
    stats::StatsType
    retcode_original::R
end

function solver_cache_init(plprob::ProfileLikelihoodProblem, target::ParameterTarget, method::OptimizationProfiler, idx, args...)

  optprob = plprob.optprob
  optpars = plprob.optpars
  #params_cache = FixedParamCache(idx, optpars[idx])
  optprob_reduced = build_optprob_reduced(optprob, optpars, idx)
  opt_cache = SciMLBase.init(optprob_reduced, get_optimizer(method); get_optimizer_opts(method)...)
  return OptimizationSolverCache(opt_cache, get_stepper(method), SciMLBase.OptimizationStats(), SciMLBase.ReturnCode.Default)
end

function solver_cache_init(plprob::ProfileLikelihoodProblem, target::FunctionTarget, method::OptimizationProfiler, idx, args...)

  #params_cache = FixedParamCache(idx, optpars[idx])
  optprob_constrained = build_optprob_constrained(plprob, idx)
  opt_cache = SciMLBase.init(optprob_constrained, get_optimizer(method); get_optimizer_opts(method)...)
  return OptimizationSolverCache(opt_cache, get_stepper(method), SciMLBase.OptimizationStats(), SciMLBase.ReturnCode.Default)
end

#get_solver_stats(solver_cache::OptimizationSolverCache) = solver_cache.opt_stats

mutable struct IntegrationSolverCache{D,O,S,R} <: AbstractSolverCache
    ode_cache::D
    opt_cache::O
    stats::S
    retcode_original::R
end

#get_solver_stats(solver_cache::IntegrationSolverCache) = solver_cache.sol.stats

function solver_cache_init(plprob::ProfileLikelihoodProblem, target::ParameterTarget, method::IntegrationProfiler, idx, dir, profile_range)
  
  ode_prob = build_odeprob_reduced(plprob, method, idx, dir, profile_range)

  # If reoptimize is requested, then create an optimization problem, create an
  # optimizer state, and register a callback.
  if method.reoptimize
    opt_method = OptimizationProfiler(; optimizer=get_optimizer(method), optimizer_opts=get_optimizer_opts(method))
    opt_solver_cache = solver_cache_init(plprob, target, opt_method, idx)
    condition(u, t, integrator) = true
    function affect!(integrator)
      set_x_fixed!(opt_solver_cache.opt_cache.reinit_cache.p, integrator.u[idx])
      opt_solver_cache.opt_cache.reinit_cache.u0 = integrator.u[1:end-1][1:end .!= idx]
      sol = solve_opt_cache(opt_solver_cache.opt_cache)
      for i in 1:length(integrator.u)-1
        i == idx && continue
        integrator.u[i] = sol[i - (i>idx)]
      end
    end
    callback = DiscreteCallback(condition, affect!)
  else
    callback = nothing
    opt_solver_cache = nothing
  end

  ode_solver_cache = SciMLBase.init(ode_prob, get_integrator(method); get_integrator_opts(method)..., callback=callback)

  return IntegrationSolverCache(ode_solver_cache, opt_solver_cache, SciMLBase.OptimizationStats(), SciMLBase.ReturnCode.Default)
end

function solver_cache_init(plprob::ProfileLikelihoodProblem, target::FunctionTarget, method::IntegrationProfiler, idx, dir, profile_range)
  
  ode_prob = build_odeprob_full(plprob, method, idx, dir, profile_range)

  # If reoptimize is requested, then create an optimization problem, create an
  # optimizer state, and register a callback.
  if method.reoptimize
    opt_method = OptimizationProfiler(; optimizer=get_optimizer(method), optimizer_opts=get_optimizer_opts(method))
    opt_solver_cache = solver_cache_init(plprob, target, opt_method, idx)
    condition(u, t, integrator) = true
    function affect!(integrator)
      opt_solver_cache.opt_cache.lcons[1] = integrator.t
      opt_solver_cache.opt_cache.ucons[1] = integrator.t
      opt_solver_cache.opt_cache.reinit_cache.u0 .= integrator.u[1:end-1]
      sol = solve_opt_cache(opt_solver_cache.opt_cache)
      for i in 1:length(integrator.u)-1
        integrator.u[i] = sol[i]
      end
    end
    callback = DiscreteCallback(condition, affect!)
  else
    callback = nothing
    opt_solver_cache = nothing
  end

  ode_solver_cache = SciMLBase.init(ode_prob, get_integrator(method); get_integrator_opts(method)..., callback=callback)

  return IntegrationSolverCache(ode_solver_cache, opt_solver_cache, SciMLBase.OptimizationStats(), SciMLBase.ReturnCode.Default)
end


mutable struct HybridSolverCache <: AbstractSolverCache

end

###################################### PROFILER CACHE ##################################

mutable struct ProfilerCache{S, Θ, SC<:AbstractSolverCache, T, ST}
  sol::S
  solver_cache::SC
  pars_cache::Θ
  pars_cur::Θ
  x_cur::T
  obj_cur::T
  profile_range::Tuple{T, T}
  iter::Int
  maxiters::Int
  verbose::Bool
  stats::ST
  retcode::Symbol
end

get_plprob(pc::ProfilerCache) = pc.sol.plprob
get_optprob(pc::ProfilerCache) = pc.sol.plprob.optprob
get_profile_target(pc::ProfilerCache) = pc.sol.plprob.target
get_profile_idx(pc::ProfilerCache) = pc.sol.idx
get_profile_dir(pc::ProfilerCache) = pc.sol.dir
get_solver_retcode(pc::ProfilerCache) = pc.solver_cache.retcode_original
get_obj_level(pc::ProfilerCache) = pc.sol.obj_level
get_pars_prev(pc::ProfilerCache) = pc.sol.pars[end-1]
get_x_prev(pc::ProfilerCache) = pc.sol.x[end-1]
get_obj_prev(pc::ProfilerCache) = pc.sol.obj[end-1]

isleft(pc::ProfilerCache) = get_profile_dir(pc) == -1
isidentifiable(pc::ProfilerCache) = pc.retcode == :Identifiable


profiler_save_values!(profiler_cache::ProfilerCache) = 
    profiler_save_values!(profiler_cache.sol, profiler_cache)

profiler_terminated(profiler_cache::ProfilerCache) = profiler_cache.retcode != :Default

function profiler_update_cache!(profiler_cache::ProfilerCache)
  @unpack x_cur, profile_range, iter, maxiters, 
  obj_cur = profiler_cache
  solver_retcode = get_solver_retcode(profiler_cache)

  dir = get_profile_dir(profiler_cache)
  obj_level = get_obj_level(profiler_cache)
  plprob = get_plprob(profiler_cache)
  if SciMLBase.successful_retcode(solver_retcode)
    profiler_save_values!(profiler_cache)
  end

  if !SciMLBase.successful_retcode(solver_retcode)
    profiler_cache.retcode = :Failure
  elseif iter > maxiters
    profiler_cache.retcode = :MaxIters
  elseif !(profile_range[1] < x_cur < profile_range[2])
    profiler_cache.retcode = :NonIdentifiable
  elseif hasthreshold(plprob) && isfinite(obj_cur) && (obj_cur ≥ obj_level)
    profiler_cache.retcode = :Identifiable
  else
    # continue profiling
    profiler_cache.retcode = :Default
  end

  return nothing
end

###################################### PROFILER VERBOSE ##################################

function init_msg(profiler_cache::ProfilerCache)
  dir = isleft(profiler_cache) ? :left : :right
  idx = get_profile_idx(profiler_cache)
  @info "Computing $dir-branch of idx=$idx profile likelihood."
end

progress_msg(profiler_cache::ProfilerCache) = progress_msg(profiler_cache, get_profile_target(profiler_cache))

function progress_msg(profiler_cache::ProfilerCache, target::ParameterTarget) 
  @info "Current parameter-$(get_profile_idx(profiler_cache)) value: $(profiler_cache.x_cur)"
end
function progress_msg(profiler_cache::ProfilerCache, target::FunctionTarget) 
  @info "Current function-$(get_profile_idx(profiler_cache)) value: $(profiler_cache.x_cur)"
end

###################################### FIXED PARAMETER CACHE ##################################

# helper type used in parameter profiling
# current implementation relies on wrapping p into FixedParamCache
# however now we can store FixedParamCache as a separate field in SolverCache
struct FixedParamCache{P}
  p::P
  idx::Base.RefValue{Int}
  x_fixed::Base.RefValue{Float64}
  gamma::Base.RefValue{Float64}
end

function FixedParamCache(p, idx::Int, x_fixed::Real, gamma::Real)
  # @assert gamma >= 0 # Todo for Sasha: justify
  _p = isnothing(p) ? SciMLBase.NullParameters() : p
  return FixedParamCache{typeof(_p)}(_p, Ref(idx), Ref(float(x_fixed)), Ref(float(gamma)))
end

get_p(p::FixedParamCache) = p.p
get_idx(p::FixedParamCache) = p.idx[]
get_x_fixed(p::FixedParamCache) = p.x_fixed[]
get_gamma(p::FixedParamCache) = p.gamma[]
set_p!(p::FixedParamCache{P}, new_p::P) where {P} = p.p = new_p
set_idx!(p::FixedParamCache, idx::Int) = p.idx[] = idx
set_x_fixed!(p::FixedParamCache, x_fixed::Float64) = p.x_fixed[] = x_fixed
set_gamma!(p::FixedParamCache, gamma::Float64) = p.gamma[] = gamma

fill_x_full!(x_full::AbstractVector{T}, x_reduced::AbstractVector{T}, p::FixedParamCache) where T<:Number = 
  fill_x_full!(x_full, x_reduced, get_idx(p), get_x_fixed(p))

function fill_x_full!(x_full::AbstractVector{T}, x_reduced::AbstractVector{T}, idx::Int, x_fixed) where T<:Number
  for i in 1:idx-1
    x_full[i] = x_reduced[i]
  end
  x_full[idx] = T(x_fixed)
  for i in idx+1:length(x_full)
    x_full[i] = x_reduced[i-1]
  end
end

function fill_x_reduced!(x_reduced::AbstractVector{T}, x_full::AbstractVector{T}, idx::Int) where T
  for i in 1:idx-1
    x_reduced[i] = x_full[i]
  end
  for i in idx+1:length(x_full)
    x_reduced[i-1] = x_full[i]
  end
end

# fill elements of reduced Matrix
function fill_x_reduced!(x_reduced::AbstractMatrix{T}, x_full::AbstractMatrix{T}, idx::Int) where T
  nr, nc = size(x_full)
  for i in 1:idx-1
    for j in 1:idx-1
      x_reduced[i, j] = x_full[i, j]
    end
    for j in idx+1:nc
      x_reduced[i, j-1] = x_full[i, j]
    end
  end

  for i in idx+1:nr
    for j in 1:idx-1
      x_reduced[i-1, j] = x_full[i, j]
    end
    for j in idx+1:nc
      x_reduced[i-1, j-1] = x_full[i, j]
    end
  end
end

