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

- `plprob::ProfileLikelihoodProblem{ParameterProfile}`: The profiling problem instance containing the parameters and likelihood function to be profiled.
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
  obj0 = evaluate_obj(_plprob.optprob, _plprob.optpars)

  !(parallel_type in (:none, :threads, :distributed)) && 
    throw(ArgumentError("Invalid `parallel_type`: $parallel_type. 
                         Supported values are :none, :threads, :distributed."))
  II = [(idx, dir) for idx in _profile_idxs(plprob) for dir in (-1, 1)]

  return __solve_parallel(_plprob, method, Val(parallel_type), II; obj0, kwargs...)
end

###################################### PARALLEL SOLVE ##################################

function __solve_parallel(plprob::ProfileLikelihoodProblem, method::AbstractProfilerMethod, ::Val{:none}, II, args...; kwargs...)

  elapsed_time = @elapsed profile_data = map(II) do I
    solve(plprob, method, I..., args...; kwargs...)
  end

  return build_profile_solution(plprob, profile_data, elapsed_time)
end

#=
function __solve_parallel(plprob::ProfileLikelihoodProblem, method::AbstractProfilerMethod, ::Val{:threads}, II; kwargs...)

  elapsed_time = @elapsed begin 
    Base.Threads.@threads for i in 1:length(II)
      idx, dir = II[i]
      profile_result = __solve_dir(plprob, method, idx, dir; kwargs...)
      output_data[i] = profile_result
    end

    profile_data = Vector{Any}(undef, length(idxs))
    for i in 1:length(idxs)
      left_profile  = output_data[2*(i-1)+1]
      right_profile = output_data[2*(i-1)+2]
      profile_data[i] = merge_profiles(left_profile, right_profile)
    end
  end
  
  return build_profile_solution(plprob, profile_data, elapsed_time)
end

function __solve_parallel(plprob::ProfileLikelihoodProblem, method::AbstractProfilerMethod, ::Val{:distributed}, II; kwargs...)

  wp = CachingPool(workers())

  elapsed_time = @elapsed begin
    output_data = Distributed.pmap(wp, II) do (idx, dir)
      solve(plprob, method, idx, dir; kwargs...)
    end

    profile_data = Vector{Any}(undef, length(idxs))
    for i in 1:length(idxs)
      left_profile  = output_data[2*(i-1)+1]
      right_profile = output_data[2*(i-1)+2]
      profile_data[i] = merge_profiles(left_profile, right_profile)
    end
  end
  
  return build_profile_solution(plprob, profile_data, elapsed_time)
end
=#

###################################### PROFILE STATE/CACHE ##################################

# mimics OptimizationState
# https://github.com/SciML/Optimization.jl/blob/master/src/state.jl
#=
"""
$(TYPEDEF)

Stores the profiler state at the current iteration 
and is passed to the callback function as the first argument.

## Fields
- `iter`: current iteration
- `x`: current solution
- `obj`: current objective value
- `grad`: current gradient
- `hess`: current hessian
"""
struct ProfilerState{X, O, G, H}
    iter::Int
    x::X
    obj::O
    grad::G
    hess::H
end

function ProfilerState(; iter = 0, x = nothing, obj = nothing,
        grad = nothing, hess = nothing)
    ProfilerState(iter, x, obj, grad, hess)
end
=#

mutable struct ProfileCurveCache
  profile_curve::C
  solver_cache::S
  idx::I
  dir::Int
  profile_bound::Float64
  iter::Int
  maxiters::Int
  verbose::Bool
  retcode::Symbol
end

# currently we use the same cache for all targets (parameters and functions), 
# but it may change to ParameterProfilerCache and FunctionProfilerCache in the future.
SciMLBase.init(plprob::ProfileLikelihoodProblem, method, idx, dir; kwargs...) = 
  SciMLBase.init(plprob, plprob.target, method, idx, dir; kwargs...)

function SciMLBase.init(plprob::ProfileLikelihoodProblem, target::AbstractProfileTarget, method::AbstractProfilerMethod, idx, dir::Int; 
  obj0 = nothing, solver_cache = nothing, maxiters = Int(1e4), verbose = false)

  x0 = evaluate_target_f(target, idx, optpars)
  if isnothing(obj0)
    verbose && @info "Computing initial values."
    _obj0 = evaluate_obj(plprob.optprob, plprob.optpars)
  else
    _obj0 = obj0
  end
  threshold = plprob.threshold
  obj_level = _obj0 + threshold

  profile_bound = dir == -1 ? _lower_bound(target, idx) : _upper_bound(target, idx)

  profile_values = ProfileCurve(Val(false), plprob, typeof(plprob.optpars),typeof(x0),typeof(_obj0), obj_level)

  if isnothing(solver_cache)
    _solver_cache = solver_init(build_scimlprob(plprob, method, idx, profile_bound), plprob, method, idx, dir, profile_bound)
  else
    # reinit_cache
    throw(ArgumentError("Providing `solver_cache` is not supported yet."))
  end

  return ProfileCurveCache{T,typeof(profile_values),typeof(method),typeof(solver_state),ReturnCode.T,typeof(idx),typeof(optpars),typeof(x0),typeof(stats)}(
    profile_values, method, solver_state, ReturnCode.Default, idx, dir, profile_bound, similar(optpars), copy(optpars), x0, obj0, 1, obj_level, maxiters, nothing, :Default, stats, verbose)


  return profile_branch_cache # ? return profiler_cache or profile_values
end


function SciMLBase.solve!(profiler_cache::AbstractProfilerCache)
  
  verbose = get_verbose(profiler_cache)
  verbose && init_msg(profiler_cache)

  profiler_save_values!(profiler_cache)
  
  while !profiler_stop(profiler_cache) 

    verbose && progress_msg(profiler_cache)
    profiler_step!(profiler_cache, get_method(profiler_cache))

    SciMLBase.successful_retcode(get_solver_retcode(profiler_cache)) && profiler_save_values!(profiler_cache)
    profiler_update_retcode!(profiler_cache)
  end
  
  profiler_finalize!(profiler_cache)
  
  return profiler_cache
end

###################################### PROFILER VERBOSE ##################################

function init_msg(profiler_cache::ProfilerState)
  dir = isleft(profiler_cache) ? :left : :right
  @info "Computing $dir-side profile"
end

function progress_msg(profiler_cache::ProfilerState{<:ParameterProfile}) 
  @info "Current parameter-$(get_idx(profiler_cache)) value: $(get_curx(profiler_cache))"
end

###################################### HELPERS ##################################

