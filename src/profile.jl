#=
function profile(plprob::PLProblem{<:FunctionProfile}, method::IntegrationProfiler, obj0::Float64, obj_level::Float64; kwargs...)
  error("Interface for profiling functions of parameters is not implemented yet...")
end
=#


"""
    profile(plprob::PLProblem, method::AbstractProfilerMethod; 
            idxs::AbstractVector{<:Int} = eachindex(get_optpars(plprob)),
            parallel_type::Symbol=:none, kwargs...)

Profiles the likelihood function for the given problem `plprob` using the specified profiling `method`.

### Arguments

- `plprob::PLProblem{ParameterProfile}`: The profiling problem instance containing the parameters and likelihood function to be profiled.
- `method::AbstractProfilerMethod`: The method to be used for profiling.
- `idxs::AbstractVector{<:Int}`: Indices of the parameters to be profiled. Defaults to all parameters. 
  (note!) Ensure that each parameter index in `idxs` has a finite `(lower, upper)` range in the PLProblem’s `profile_range` – the profile procedure will validate this.
- `parallel_type::Symbol`: Specifies the type of parallelism to be used. Supported values: `:none, :threads, :distributed`. Defaults to `:none`.
- `maxiters::Int`: Maximum number of iterations for one branch (left and right) of the profiling process. Defaults to `1e4`.
- `verbose::Bool`: Indicates whether to display the progress of the profiling process. Defaults to `false`.

### Returns

- Returns the profiling results `PLSolution`.

### Example

```julia
plprob = PLProblem(optprob, optpars, [(-10.,10.), (-5.,5.)])
method = OptimizationProfiler(optimizer = Optimization.LBFGS(), stepper = FixedStep())
sol = profile(plprob, method; idxs=[1])
```
"""
function profile(plprob::PLProblem{ParameterProfile}, method::AbstractProfilerMethod; 
  idxs::AbstractVector{<:Int} = eachindex(get_optpars(plprob)),
  parallel_type::Symbol=:none, kwargs...)

  # validation
  @assert parallel_type in (:none, :threads, :distributed)
  optpars = get_optpars(plprob)
  checkbounds(optpars, idxs)
  validate_profile_range(plprob, idxs)

  return __profile(plprob, method, Val(parallel_type), idxs; kwargs...)
end


function __profile(plprob::PLProblem, method::AbstractProfilerMethod, ::Val{:none}, idxs; kwargs...)

  elapsed_time = @elapsed profile_data = map(idxs) do idx
    left_profile  = __profile_dir(plprob, method, idx, -1; kwargs...)
    right_profile = __profile_dir(plprob, method, idx,  1; kwargs...)
    merge_profiles(left_profile, right_profile)
  end

  return build_profile_solution(plprob, profile_data, elapsed_time)
end

function __profile(plprob::PLProblem, method::AbstractProfilerMethod, ::Val{:threads}, idxs; kwargs...)

  input_data = [(idx, dir) for idx in idxs for dir in (-1, 1)]
  output_data = Vector{Any}(undef, length(input_data))

  elapsed_time = @elapsed begin 
    Base.Threads.@threads for i in 1:length(input_data)
      idx, dir = input_data[i]
      profile_result = __profile_dir(plprob, method, idx, dir; kwargs...)
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

function __profile(plprob::PLProblem, method::AbstractProfilerMethod, ::Val{:distributed}, idxs; kwargs...)

  input_data = [(idx, dir) for idx in idxs for dir in (-1, 1)]

  elapsed_time = @elapsed begin
    output_data = Distributed.pmap(input_data) do (idx, dir)
      __profile_dir(plprob, method, idx, dir; kwargs...)
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

function __profile_dir(plprob::PLProblem, method::AbstractProfilerMethod, idx::Int, dir::Int; kwargs...)
  
  profiler_state = profiler_init(plprob, method, idx, dir; kwargs...)
  __profile_dir!(profiler_state)

  return profiler_state # ? return profiler_state or profile_values
end


function __profile_dir!(profiler_state::ProfilerState)
  
  verbose = get_verbose(profiler_state)
  verbose && init_msg(profiler_state)

  profiler_save_values!(profiler_state)
  
  while !profiler_stop(profiler_state) 

    verbose && progress_msg(profiler_state)
    profiler_step!(profiler_state, get_method(profiler_state))

    SciMLBase.successful_retcode(get_solver_retcode(profiler_state)) && profiler_save_values!(profiler_state)
    profiler_update_retcode!(profiler_state)
  end
  
  profiler_finalize!(profiler_state)
  
  return profiler_state
end

###################################### PROFILER VERBOSE ##################################

function init_msg(profiler_state::ProfilerState)
  dir = isleft(profiler_state) ? :left : :right
  @info "Computing $dir-side profile"
end

function progress_msg(profiler_state::ProfilerState{<:ParameterProfile}) 
  @info "Current parameter-$(get_idx(profiler_state)) value: $(get_curx(profiler_state))"
end

###################################### HELPERS ##################################

function validate_profile_range(plprob::PLProblem{<:ParameterProfile}, idxs::AbstractVector{<:Int})
  optpars = get_optpars(plprob)
  profile_range = get_profile_range(plprob)
  
  if profile_range isa Tuple
    for idx in idxs
      validate_profile_range(optpars[idx], profile_range)
    end
  else
    for idx in idxs
      validate_profile_range(optpars[idx], profile_range[idx]) 
    end
  end
end

function validate_profile_range(x::Number, profile_range)
  !(profile_range isa Tuple) &&
    throw(ArgumentError("`profile_range` must contain finite `(lower, upper)` tuples for all parameters selected for profiling."))
  lb, ub = profile_range
  validate_profile_bound(lb)
  validate_profile_bound(ub)
  !(lb <= x <= ub) &&
    throw(ArgumentError("The initial values provided for profiling parameters must lie within the specified `profile_range`: `profile_range[idx][1] ≤ x[idx] ≤ profile_range[idx][2]`"))
end

function validate_profile_bound(bound) 
  (isnothing(bound) || isinf(bound)) &&
    throw(ArgumentError("Each parameter selected for profiling must have a finite `(lower, upper)` range specified in `profile_range`"))
end
