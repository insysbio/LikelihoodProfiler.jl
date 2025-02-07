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
- `parallel_type::Symbol`: Specifies the type of parallelism to be used. Defaults to `:none`.
- `kwargs...`: Additional keyword arguments to be passed to the profiling method.

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
   
  optpars = get_optpars(plprob)
  checkbounds(optpars, idxs)

  return __profile(plprob, method, Val(parallel_type), idxs; kwargs...)
end


function __profile(plprob::PLProblem, method::AbstractProfilerMethod, ::Val{:none}, idxs; kwargs...)

  sciml_prob = build_scimlprob(plprob, method)
  #solver_state = solver_init(plprob, method)

  elapsed_time = @elapsed profile_data = map(idxs) do idx
    left_profile  = __profile_dir(plprob, method, sciml_prob, idx, -1; kwargs...)
    right_profile = __profile_dir(plprob, method, sciml_prob, idx,  1; kwargs...)
    merge_profiles(left_profile, right_profile)
  end

  return build_profile_solution(plprob, profile_data, elapsed_time)
end



function __profile_dir(plprob::PLProblem, method::AbstractProfilerMethod, sciml_prob, idx::Int, dir::Int; kwargs...)
  
  profiler_state = profiler_init(plprob, method, sciml_prob, idx, dir; kwargs...)
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
#=
function progress_msg(profiler_state::ProfilerState{<:FunctionProfile})
  @info "Current profile function value: $(get_curx(profiler_state))"
end
=#