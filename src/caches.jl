mutable struct OptimizationProfilerCache{O,StepType,StatsType,R} <: AbstractProfilerMethodCache
    opt_cache::O
    stepper::StepType
    stats::StatsType
    retcode_original::R
end

function init_method_cache(plprob::ProfileLikelihoodProblem, target::ParameterTarget, method::OptimizationProfiler, idx, args...)

  optprob = plprob.optprob
  optpars = plprob.optpars
  #params_cache = FixedParamCache(idx, optpars[idx])
  optprob_reduced = build_optprob_reduced(optprob, optpars, idx)
  opt_cache = SciMLBase.init(optprob_reduced, get_optimizer(method); get_optimizer_opts(method)...)
  return OptimizationProfilerCache(opt_cache, get_stepper(method), SciMLBase.OptimizationStats(), SciMLBase.ReturnCode.Default)
end

function init_method_cache(plprob::ProfileLikelihoodProblem, target::FunctionTarget, method::OptimizationProfiler, idx, args...)

  #params_cache = FixedParamCache(idx, optpars[idx])
  optprob_constrained = build_optprob_constrained(plprob, idx)
  opt_cache = SciMLBase.init(optprob_constrained, get_optimizer(method); get_optimizer_opts(method)...)
  return OptimizationProfilerCache(opt_cache, get_stepper(method), SciMLBase.OptimizationStats(), SciMLBase.ReturnCode.Default)
end

#get_solver_stats(method_cache::OptimizationProfilerCache) = method_cache.opt_stats

mutable struct IntegrationProfilerCache{D,O,S,R} <: AbstractProfilerMethodCache
    ode_cache::D
    opt_cache::O
    stats::S
    retcode_original::R
end

#get_solver_stats(method_cache::IntegrationProfilerCache) = method_cache.sol.stats

function init_method_cache(plprob::ProfileLikelihoodProblem, target::ParameterTarget, method::IntegrationProfiler, idx, dir, profile_range)
  
  ode_prob = build_odeprob_reduced(plprob, method, idx, dir, profile_range)

  # If reoptimize is requested, then create an optimization problem, create an
  # optimizer state, and register a callback.
  if method.reoptimize
    opt_method = OptimizationProfiler(; optimizer=get_optimizer(method), optimizer_opts=get_optimizer_opts(method))
    opt_method_cache = init_method_cache(plprob, target, opt_method, idx)
    condition(u, t, integrator) = true
    function affect!(integrator)
      set_x_fixed!(opt_method_cache.opt_cache.reinit_cache.p, integrator.u[idx])
      fill_x_reduced!(opt_method_cache.opt_cache.reinit_cache.u0, integrator.u.theta, idx)
      sol = solve_opt_cache(opt_method_cache.opt_cache)

      if !SciMLBase.successful_retcode(sol.retcode)
        @warn "Reoptimization returned $(sol.retcode) at profile point x = $(integrator.t). Profiling is interrupted."
        SciMLBase.terminate!(integrator, sol.retcode)
        return nothing
      end

      for i in 1:length(integrator.u)-1
        i == idx && continue
        integrator.u[i] = sol[i - (i > idx)]
      end
      return nothing
    end
    callback = SciMLBase.DiscreteCallback(condition, affect!)
  else
    callback = nothing
    opt_method_cache = nothing
  end

  ode_cache = SciMLBase.init(ode_prob, get_integrator(method); get_integrator_opts(method)..., callback=callback)

  return IntegrationProfilerCache(ode_cache, opt_method_cache, SciMLBase.OptimizationStats(), SciMLBase.ReturnCode.Default)
end

function init_method_cache(plprob::ProfileLikelihoodProblem, target::FunctionTarget, method::IntegrationProfiler, idx, dir, profile_range)
  
  ode_prob = build_odeprob_full(plprob, method, idx, dir, profile_range)

  # If reoptimize is requested, then create an optimization problem, create an
  # optimizer state, and register a callback.
  if method.reoptimize
    opt_method = OptimizationProfiler(; optimizer=get_optimizer(method), optimizer_opts=get_optimizer_opts(method))
    opt_method_cache = init_method_cache(plprob, target, opt_method, idx)
    condition(u, t, integrator) = true
    function affect!(integrator)
      opt_method_cache.opt_cache.lcons[1] = integrator.t
      opt_method_cache.opt_cache.ucons[1] = integrator.t
      copyto!(opt_method_cache.opt_cache.reinit_cache.u0, integrator.u.theta)
      sol = solve_opt_cache(opt_method_cache.opt_cache)

      if !SciMLBase.successful_retcode(sol.retcode)
        @warn "Reoptimization returned $(sol.retcode) at profile point x = $(integrator.t). Profiling is interrupted."
        SciMLBase.terminate!(integrator, sol.retcode)
        return nothing
      end

      for i in 1:length(integrator.u)-1
        integrator.u[i] = sol[i]
      end
      if get_matrix_type(method) == :hessian
        update_lagrange_multiplier!(integrator.u, plprob.optprob, plprob.target.fs[idx])
      end
      return nothing
    end
    callback = SciMLBase.DiscreteCallback(condition, affect!)
  else
    callback = nothing
    opt_method_cache = nothing
  end

  ode_cache = SciMLBase.init(ode_prob, get_integrator(method); get_integrator_opts(method)..., callback=callback)

  return IntegrationProfilerCache(ode_cache, opt_method_cache, SciMLBase.OptimizationStats(), SciMLBase.ReturnCode.Default)
end


mutable struct HybridProfilerCache <: AbstractProfilerMethodCache

end

###################################### PROFILE BRANCH CACHE ##################################

mutable struct ProfileBranchCache{S, Θ, MC<:AbstractProfilerMethodCache, T, ST}
  sol::S
  method_cache::MC
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

get_plprob(pc::ProfileBranchCache) = pc.sol.plprob
get_threshold(pc::ProfileBranchCache) = pc.sol.plprob.threshold
get_optpars(pc::ProfileBranchCache) = pc.sol.plprob.optpars
get_optprob(pc::ProfileBranchCache) = pc.sol.plprob.optprob
get_profile_target(pc::ProfileBranchCache) = pc.sol.plprob.target
get_profile_idx(pc::ProfileBranchCache) = pc.sol.idx
get_profile_dir(pc::ProfileBranchCache) = pc.sol.dir
get_solver_retcode(pc::ProfileBranchCache) = pc.method_cache.retcode_original
get_obj_level(pc::ProfileBranchCache) = pc.sol.obj_level
get_pars_prev(pc::ProfileBranchCache) = pc.sol.pars[end-1]
get_x_prev(pc::ProfileBranchCache) = pc.sol.x[end-1]
get_obj_prev(pc::ProfileBranchCache) = pc.sol.obj[end-1]

isleft(pc::ProfileBranchCache) = get_profile_dir(pc) == -1
isidentifiable(pc::ProfileBranchCache) = pc.retcode == :Identifiable

evaluate_obj(cache::ProfileBranchCache, u=cache.pars_cur) =
    evaluate_obj(get_plprob(cache), u)
    
function evaluate_target_f(cache::ProfileBranchCache, u=cache.pars_cur)
    plprob = get_plprob(cache)
    return evaluate_target_f(
        plprob.target, get_profile_idx(cache), u, plprob.optprob.p
    )
end

profiler_save_values!(branch_cache::ProfileBranchCache) =
    profiler_save_values!(branch_cache.sol, branch_cache)

profiler_terminated(branch_cache::ProfileBranchCache) = branch_cache.retcode != :Default

function profiler_update_cache!(branch_cache::ProfileBranchCache)
  @unpack x_cur, profile_range, iter, maxiters, 
  obj_cur = branch_cache
  solver_retcode = get_solver_retcode(branch_cache)
  
  plprob = get_plprob(branch_cache)
  if SciMLBase.successful_retcode(solver_retcode)
    profiler_save_values!(branch_cache)
  end

  if !SciMLBase.successful_retcode(solver_retcode)
    branch_cache.retcode = :Failure
  elseif iter > maxiters
    branch_cache.retcode = :MaxIters
  elseif hasthreshold(plprob) && isfinite(obj_cur) && reached_threshold(branch_cache)
    branch_cache.retcode = :Identifiable
  elseif reached_profile_bound(branch_cache)
    branch_cache.retcode = :NonIdentifiable
  else
    # continue profiling
    branch_cache.retcode = :Default
  end

  return nothing
end

reached_threshold(branch_cache::ProfileBranchCache) = reached_threshold(branch_cache, branch_cache.method_cache)

function reached_threshold(branch_cache::ProfileBranchCache, ::AbstractProfilerMethodCache)
   @unpack obj_cur = branch_cache
   obj_level = get_obj_level(branch_cache)
   return obj_cur ≥ obj_level
end

reached_profile_bound(pc::ProfileBranchCache) = reached_profile_bound(pc, pc.method_cache)

function reached_profile_bound(pc::ProfileBranchCache, ::AbstractProfilerMethodCache)
    return !(pc.profile_range[1] < pc.x_cur < pc.profile_range[2])
end

###################################### PROFILER VERBOSE ##################################

function init_msg(branch_cache::ProfileBranchCache)
  dir = isleft(branch_cache) ? :left : :right
  idx = get_profile_idx(branch_cache)
  @info "Computing $dir-branch of idx=$idx profile likelihood."
end

progress_msg(branch_cache::ProfileBranchCache) = progress_msg(branch_cache, get_profile_target(branch_cache))

function progress_msg(branch_cache::ProfileBranchCache, target::ParameterTarget)
  @info "Current parameter-$(get_profile_idx(branch_cache)) value: $(branch_cache.x_cur)"
end
function progress_msg(branch_cache::ProfileBranchCache, target::FunctionTarget)
  @info "Current function-$(get_profile_idx(branch_cache)) value: $(branch_cache.x_cur)"
end

###################################### FIXED PARAMETER CACHE ##################################

# helper type used in parameter profiling
# current implementation relies on wrapping p into FixedParamCache
# however now we can store FixedParamCache as a separate field in the method cache
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

