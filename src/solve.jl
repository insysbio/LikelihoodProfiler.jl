#=
function SciMLBase.solve(plprob::ProfileLikelihoodProblem{<:FunctionProfile}, method::IntegrationProfiler, obj0::Float64, obj_level::Float64; kwargs...)
  error("Interface for profiling functions of parameters is not implemented yet...")
end
=#


"""
    solve(plprob::ProfileLikelihoodProblem, method::AbstractProfilerMethod; 
            parallel_type::Symbol=:none, maxiters::Int=1e4, verbose::Bool=false)

Profiles the likelihood function for the given problem `plprob` using the specified profiling `method`.

### Arguments

- `plprob::ProfileLikelihoodProblem`: The profiling problem instance containing the parameters and likelihood function to be profiled.
- `method::AbstractProfilerMethod`: The method to be used for profiling.
- `reoptimize_init::Bool=false`: If `true`, re-optimizes the model at the provided initial parameter values `optpars` before profiling. Defaults to `false`.
- `parallel_type::Symbol`: Specifies the type of parallelism to be used. Supported values: `:none, :threads, :distributed`. Defaults to `:none`.
- `maxiters::Int`: Maximum number of iterations for one branch (left and right) of the profiling process. Defaults to `1e4`.
- `verbose::Bool`: Indicates whether to display the progress of the profiling process. Defaults to `false`.

### Returns

- Returns the profiling results `ProfileLikelihoodSolution`.

### Example

```julia
plprob = ProfileLikelihoodProblem(optprob, optpars; idxs=1, profile_lower=-10., profile_upper=10.)
method = OptimizationProfiler(optimizer = Optimization.LBFGS(), stepper = FixedStep())
sol = solve(plprob, method)
```
"""
function SciMLBase.solve(plprob::ProfileLikelihoodProblem, method::AbstractProfilerMethod; 
  reoptimize_init::Bool=false, parallel_type::Symbol=:none, verbose::Bool=false, kwargs...)

  check_prob_alg(plprob, method)

  if reoptimize_init
    !has_optimizer(method) && 
      throw(ArgumentError("`method` must have a valid `optimizer` provided when `reoptimize_init=true`."))

    # start from user 'optpars'
    optprob0 = remake(plprob.optprob; u0 = plprob.optpars)
    s = solve(optprob0, _optimizer(method); _optimizer_opts(method)...)
    if !SciMLBase.successful_retcode(s)
      @warn "Re-optimization at initial parameter values returned $(s.retcode). Proceeding with the provided initial parameters."
      _plprob = plprob
    else
      _plprob = remake(plprob; optpars=s.u)
    end
  else
    _plprob = plprob
  end

  verbose && @info "Computing initial values."
  obj0 = evaluate_obj(_plprob, _plprob.optpars)

  !(parallel_type in (:none, :threads, :distributed)) && 
    throw(ArgumentError("Invalid `parallel_type`: $parallel_type. 
                         Supported values are :none, :threads, :distributed."))
  II = [(idx, dir) for idx in get_profile_idxs(plprob.target) for dir in (-1, 1)]

  return __solve_parallel(_plprob, method, Val(parallel_type), II; obj0, verbose, kwargs...)
end

###################################### PARALLEL SOLVE ##################################

function __solve_parallel(plprob::ProfileLikelihoodProblem, method::AbstractProfilerMethod, ::Val{:none}, II, args...; kwargs...)

  elapsed_time = @elapsed profile_data = map(II) do I
    solve(plprob, method, I..., args...; kwargs...)
  end

  return build_profile_solution(plprob, profile_data, elapsed_time)
end


function __solve_parallel(plprob::ProfileLikelihoodProblem, method::AbstractProfilerMethod, ::Val{:threads}, II, args...; kwargs...)

  elapsed_time = @elapsed begin 
    profile_data = Vector{Any}(undef, length(II))
    Base.Threads.@threads for i in 1:length(II)
      profile_data[i] = solve(plprob, method, II[i]..., args...; kwargs...)
    end
  end
  
  return build_profile_solution(plprob, profile_data, elapsed_time)
end

function __solve_parallel(plprob::ProfileLikelihoodProblem, method::AbstractProfilerMethod, ::Val{:distributed}, II, args...; kwargs...)

  wp = CachingPool(workers())

  elapsed_time = @elapsed profile_data = Distributed.pmap(wp, II) do I
    solve(plprob, method, I..., args...; kwargs...)
  end
  
  return build_profile_solution(plprob, profile_data, elapsed_time)
end


###################################### PROFILE STATE/CACHE ##################################

# currently we use the same cache for all targets (parameters and functions), 
# but it may change to ParameterProfilerCache and FunctionProfilerCache in the future.
SciMLBase.init(plprob::ProfileLikelihoodProblem, method, idx, dir; kwargs...) = 
  SciMLBase.init(plprob, plprob.target, method, idx, dir; kwargs...)

function SciMLBase.init(plprob::ProfileLikelihoodProblem, target::AbstractProfileTarget, method::AbstractProfilerMethod, idx, dir::Int; 
  obj0 = nothing, solver_cache = nothing, maxiters = Int(1e4), verbose = false)

  (dir == -1 || dir == +1) ||
        throw(ArgumentError("Profile direction `dir` must be -1 or +1"))
  #=
  n = length(plprob.optpars)
    (1 <= idx <= n) ||
        throw(ArgumentError("`idx` must be within 1:$n"))
  =#
  θ0 = plprob.optpars
  x0 = evaluate_target_f(target, idx, θ0)
  _obj0 = isnothing(obj0) ? evaluate_obj(plprob.optprob, θ0) : obj0

  τ  = plprob.threshold
  obj_level = _obj0 + τ

  lb = get_profile_lb(target, idx); ub = get_profile_ub(target, idx)

  T = float(promote_type(eltype(θ0), typeof(x0), typeof(_obj0), typeof(τ), typeof(lb), typeof(ub)))
  θ0_typed = T.(θ0)
  profile_range = (T(lb), T(ub))
  sol = solution_init(plprob, idx, dir, θ0_typed, T(x0), T(_obj0), T(obj_level))

  if isnothing(solver_cache)
    _solver_cache = solver_cache_init(plprob, method, idx, dir, profile_range)
  else
    # reinit_cache
    throw(ArgumentError("Providing `solver_cache` is not supported yet."))
  end

  return profiler_init(sol, _solver_cache, θ0_typed, T(x0), T(_obj0), profile_range, Int(maxiters), verbose)
end


function SciMLBase.solve!(profiler_cache::ProfilerCache)

  verbose = profiler_cache.verbose
  verbose && init_msg(profiler_cache)
  while !profiler_terminated(profiler_cache) 

    verbose && progress_msg(profiler_cache)
    profiler_step!(profiler_cache)

    profiler_update_cache!(profiler_cache)
  end
  
  profiler_finalize_solution!(profiler_cache)
  
  return profiler_cache.sol
end

solver_cache_init(plprob, method, idx, dir, profile_range) = solver_cache_init(plprob, plprob.target, method, idx, dir, profile_range)
profiler_step!(profiler_cache::ProfilerCache) = profiler_step!(profiler_cache, get_profile_target(profiler_cache), profiler_cache.solver_cache)
profiler_finalize_solution!(profiler_cache::ProfilerCache) = profiler_finalize_solution!(profiler_cache::ProfilerCache, profiler_cache.sol)


################################### CHECK PROB ALG ##################################

const OPTF_GRAD_ERROR = "The algorithm you are using requires gradient, 
  but the provided `OptimizationFunction` has neither a user-supplied gradient nor an AD tag.
  Provide a gradient (OptimizationFunction(...; grad=...)) or set an AD tag, e.g. `adtype=AutoForwardDiff()`."

const OPTF_HESS_ERROR = "The algorithm you are using requires Hessian,
  but the provided `OptimizationFunction` has neither a user-supplied Hessian nor a Hessian-capable AD tag.
  Provide a Hessian `(OptimizationFunction(...; hess=...))` or choose an AD that can compute Hessians, e.g. `AutoForwardDiff()` or `AutoFiniteDiff()`."

check_prob_alg(plprob::ProfileLikelihoodProblem, method::AbstractProfilerMethod) = 
  check_prob_alg(plprob, plprob.target, method)

check_prob_alg(plprob::ProfileLikelihoodProblem, target::AbstractProfileTarget, method::AbstractProfilerMethod) = nothing

check_prob_alg(plprob::ProfileLikelihoodProblem, target::FunctionTarget, method::OptimizationProfiler) = 
  throw(ArgumentError("`OptimizationProfiler` currently doesn't support function profiling. Please consider using other profiling algorithms."))

function check_prob_alg(plprob::ProfileLikelihoodProblem, target::AbstractProfileTarget, method::OptimizationProfiler)
  optf = plprob.optprob.f
  opt_alg = method.optimizer

  SciMLBase.requiresgradient(opt_alg) && !hasgrad(optf) &&
    throw(ArgumentError(OPTF_GRAD_ERROR))

  SciMLBase.requireshessian(opt_alg) && !hashess(optf) &&
    throw(ArgumentError(OPTF_HESS_ERROR))
  return nothing
end

function check_prob_alg(plprob::ProfileLikelihoodProblem, target::ParameterTarget, method::IntegrationProfiler)
  optf = plprob.optprob.f

  if method.matrix_type == :hessian
    !hashess(optf) && throw(ArgumentError(OPTF_HESS_ERROR))
  else
    !hasgrad(optf) && throw(ArgumentError(OPTF_GRAD_ERROR))
  end

  if method.reoptimize
    opt_alg = method.optimizer
    SciMLBase.requiresgradient(opt_alg) && !hasgrad(optf) &&
      throw(ArgumentError(OPTF_GRAD_ERROR))

    SciMLBase.requireshessian(opt_alg) && !hashess(optf) &&
      throw(ArgumentError(OPTF_HESS_ERROR))
  end
  return nothing
end

function check_prob_alg(plprob::ProfileLikelihoodProblem, target::FunctionTarget, method::IntegrationProfiler)
  optf = plprob.optprob.f
  profile_fs = get_profile_fs(target)

  if method.matrix_type == :hessian
    !hashess(optf) && throw(ArgumentError(OPTF_HESS_ERROR))
    for f in profile_fs
      !hashess(f) && throw(ArgumentError(OPTF_HESS_ERROR))
      !hasgrad(f) && throw(ArgumentError(OPTF_GRAD_ERROR))
    end
  else
    !hasgrad(optf) && throw(ArgumentError(OPTF_GRAD_ERROR))
    for f in profile_fs
      !hasgrad(f) && throw(ArgumentError(OPTF_GRAD_ERROR))
    end
  end

  if method.reoptimize
    opt_alg = method.optimizer
    !SciMLBase.allowsconstraints(opt_alg) && 
      throw(ArgumentError("`IntegrationProfiler` with `reoptimize=true` requires an optimization algorithm that supports constraints."))
    SciMLBase.requiresgradient(opt_alg) && !hasgrad(optf) &&
      throw(ArgumentError(OPTF_GRAD_ERROR))

    SciMLBase.requireshessian(opt_alg) && !hashess(optf) &&
      throw(ArgumentError(OPTF_HESS_ERROR))
  end
  return nothing
end

hasgrad(f::OptimizationFunction) = 
  (hasfield(typeof(f), :grad)  && f.grad  !== nothing) ||
  (hasfield(typeof(f), :adtype) && !(f.adtype isa SciMLBase.NoAD))

hashess(f::OptimizationFunction) = 
  (hasfield(typeof(f), :hess)  && f.hess  !== nothing) ||
  (hasfield(typeof(f), :adtype) && !(f.adtype isa SciMLBase.NoAD))

