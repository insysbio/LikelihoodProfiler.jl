
mutable struct ProfilerState{T,pvType,methodType,SS,SR,ID,parType,xType,statsType}
  profile_values::pvType
  method::methodType
  solver_state::SS
  solver_retcode::SR
  idx::ID
  dir::Int
  profile_bound::Float64
  parscache::parType
  curpars::parType
  curx::xType
  curobj::Float64
  numiter::Int
  #curgrad::gradType
  #curhess::hessType
  obj_level::Float64
  maxiters::Int
  endpoint::Union{Nothing,Float64}
  retcode::Symbol
  stats::statsType
  verbose::Bool
end

isleft(ps::ProfilerState) = get_dir(ps::ProfilerState) == -1
isidentifiable(ps::ProfilerState) = get_retcode(ps::ProfilerState) == :Identifiable
############################### SELECTORS ###############################

get_profile_values(ps::ProfilerState) = ps.profile_values
get_plprob(ps::ProfilerState) = get_plprob(get_profile_values(ps))
get_optprob(ps::ProfilerState) = get_optprob(get_plprob(ps))
get_method(ps::ProfilerState) = ps.method
get_solver_state(ps::ProfilerState) = ps.solver_state
get_solver_retcode(ps::ProfilerState) = ps.solver_retcode
#get_plprob(ps::ProfilerState) = ps.plprob
#get_parfunc(ps::ProfilerState) = get_parfunc(get_plprob(ps))
get_idx(ps::ProfilerState) = ps.idx
get_dir(ps::ProfilerState) = ps.dir
get_profile_bound(ps::ProfilerState) = ps.profile_bound
get_parscache(ps::ProfilerState) = ps.parscache
get_curpars(ps::ProfilerState) = ps.curpars
get_curx(ps::ProfilerState) = ps.curx
get_curobj(ps::ProfilerState) = ps.curobj
get_numiter(ps::ProfilerState) = ps.numiter
# get_curgrad(ps::ProfilerState) = ps.curgrad
# get_curhess(ps::ProfilerState) = ps.curhess
get_obj_level(ps::ProfilerState) = ps.obj_level
get_maxiters(ps::ProfilerState) = ps.maxiters
get_endpoint(ps::ProfilerState) = ps.endpoint
get_retcode(ps::ProfilerState) = ps.retcode
get_stats(ps::ProfilerState) = ps.stats
get_verbose(ps::ProfilerState) = ps.verbose

get_profile_type(ps::ProfilerState{T}) where T= T
#get_optprob(ps::ProfilerState) = get_optprob(get_plprob(ps))
#get_optpars(ps::ProfilerState) = get_optpars(get_plprob(ps))
#get_optobj(ps::ProfilerState) = get_optobj(get_plprob(ps))
#get_profile_range(ps::ProfilerState) = get_profile_range(get_plprob(ps))
#get_obj_level(ps::ProfilerState) = get_obj_level(get_plprob(ps))


profiler_stop(profiler::ProfilerState) = get_retcode(profiler) != :Default

function profiler_update_retcode!(profiler_state::ProfilerState)

  dir = get_dir(profiler_state)
  curx = get_curx(profiler_state)
  profile_bound = get_profile_bound(profiler_state)
  numiter = get_numiter(profiler_state)
  maxiters = get_maxiters(profiler_state)
  curobj = get_curobj(profiler_state)
  obj_level = get_obj_level(profiler_state)
  solver_retcode = get_solver_retcode(profiler_state)
  
  if !SciMLBase.successful_retcode(solver_retcode)
    profiler_state.retcode = :Failed
  elseif numiter > maxiters
    profiler_state.retcode = :MaxIters
  elseif dir*curx >= dir*profile_bound
    profiler_state.retcode = :NonIdentifiable
  elseif curobj >= obj_level
    profiler_state.retcode = :Identifiable
  else
    # continue profiling
    profiler_state.retcode = :Default
  end
end


function profiler_finalize!(profiler_state::ProfilerState)
  profile_values = get_profile_values(profiler_state)
  isidentifiable(profiler_state) && (profiler_state.endpoint = interpolate_endpoint(profile_values))
  isleft(profiler_state) && reverse_profile_values!(profile_values)
  return nothing
end

#=
function profiler_finalize!(profiler_state::ProfilerState)
  profile_values = get_profile_values(profiler_state)
  ep = isidentifiable(profiler_state) ? interpolate_endpoint(profile_values) : nothing
  if isleft(profiler_state)
    endpoints = (lower=ep,)
    retcodes = (lower=,)
    stats = (lower=ep,)
    reverse_profile_values!(profile_values)
  else
    endpoints = (upper=ep,)
    retcodes = (upper=,)
    stats = (lower=ep,)
  end
  return nothing
end
=#

function profiler_save_values!(profiler_state::ProfilerState)
  profile_values = get_profile_values(profiler_state)
  numiter = get_numiter(profiler_state)
  SciMLBase.copyat_or_push!(profile_values.pars, numiter, get_curpars(profiler_state))
  SciMLBase.copyat_or_push!(profile_values.x,    numiter, get_curx(profiler_state))
  SciMLBase.copyat_or_push!(profile_values.obj,  numiter, get_curobj(profiler_state))
  return nothing
end

#reverse_profile_values!(ps::ProfilerState) = reverse_profile_values!(get_profile_values(ps))
function reverse_profile_values!(pv::ProfileValues)
  reverse!(pv.pars)
  reverse!(pv.x)
  reverse!(pv.obj)
  return nothing
end

function profiler_init(plprob::PLProblem{T}, method::AbstractProfilerMethod, sciml_prob, idx, dir;
  maxiters = Int(1e4), verbose = false) where T
  
  verbose && @info "Computing initial values."
  optprob = get_optprob(plprob)
  optpars = get_optpars(plprob)
  threshold = get_threshold(plprob)
  
  obj0 = compute_optf(optprob, optpars)
  obj_level = obj0 + threshold
  x0 = optpars[idx]

  profile_range = get_profile_range(plprob)
  profile_lb, profile_ub = profile_range isa Tuple ? profile_range : profile_range[idx]
  profile_bound = dir == -1 ? profile_lb : profile_ub

  solver_state = solver_init(sciml_prob, plprob, method, idx, dir, profile_bound)
  #solver_reinit!(solver_state, plprob, method, idx, dir, profile_bound)
  profile_values = ProfileValues(Val(false), plprob, typeof(optpars),typeof(x0),typeof(obj0), obj_level)

  stats = get_solver_stats(solver_state)
  return ProfilerState{T,typeof(profile_values),typeof(method),typeof(solver_state),ReturnCode.T,typeof(idx),typeof(optpars),typeof(x0),typeof(stats)}(
    profile_values, method, solver_state, ReturnCode.Default, idx, dir, profile_bound, similar(optpars), copy(optpars), x0, obj0, 1, obj_level, maxiters, nothing, :Default, stats, verbose)
end

#=
function profiler_init(plprob::PLProblem{T}, method::IntegrationProfiler, solver_state, idx, dir, profile_bound, obj0, obj_level; maxiters = Int(1e4)) where T
  
  optpars = get_optpars(plprob)
  x0 = optpars[idx]

  #length(optpars)>1 && 
  solver_state = solver_reinit!(solver_state, plprob, method, idx, dir, profile_bound)
  profile_values = ProfileValues(Val(false), plprob, typeof(optpars),typeof(x0),typeof(obj0), obj_level)

  stats = get_solver_stats(solver_state)
  return ProfilerState{T,typeof(profile_values),typeof(method),typeof(solver_state),ReturnCode.T,typeof(idx),typeof(optpars),typeof(x0),typeof(stats)}(
    profile_values, method, solver_state, ReturnCode.Default, idx, dir, profile_bound, similar(optpars), copy(optpars), x0, obj0, 1, obj_level, maxiters, nothing, :Default, stats)
end
=#
function merge_profiles(left_profile::ProfilerState, right_profile::ProfilerState)
  left_profile_values = get_profile_values(left_profile)
  right_profile_values = get_profile_values(right_profile)
  plprob = get_plprob(left_profile_values)
  obj_level = get_obj_level(left_profile_values)

  pars = [left_profile_values.pars; right_profile_values.pars]
  x = [left_profile_values.x; right_profile_values.x]
  obj = [left_profile_values.obj; right_profile_values.obj]
  retcodes = (left_profile.retcode, right_profile.retcode)
  endpoints = (left_profile.endpoint, right_profile.endpoint)
  stats = (left_profile.stats, right_profile.stats)

  return ProfileValues(Val(false), plprob, pars, x, obj, obj_level, retcodes, endpoints, stats)
end

function merge_profiles(left_profile::ProfileValues, right_profile::ProfileValues)
  plprob = get_plprob(left_profile)
  obj_level = get_obj_level(left_profile)

  pars = [left_profile.pars; right_profile.pars]
  x = [left_profile.x; right_profile.x]
  obj = [left_profile.obj; right_profile.obj]
  retcodes = (left_profile.retcodes[1], right_profile.retcodes[2])
  endpoints = (left_profile.endpoints[1], right_profile.endpoints[2])
  stats = (left_profile.stats[1], right_profile.stats[2])

  return ProfileValues(Val(false), plprob, pars, x, obj, obj_level, retcodes, endpoints, stats)
end

function profiler_step!(profiler::ProfilerState, method::OptimizationProfiler)
  stepper = get_stepper(method)

  pars_guess = compute_next_pars!(profiler, stepper)
  sol = solve_optcache(profiler)

  if SciMLBase.successful_retcode(sol.retcode)
    idx = get_idx(profiler)
    fill_x_full!(profiler.curpars, sol.u, idx, pars_guess[idx])
    profiler.curx = pars_guess[idx]
    profiler.curobj = sol.objective
    profiler.solver_retcode = sol.retcode
    profiler.numiter += 1
  else
    @warn "Solver returned $(sol.retcode) retcode at profile point x = $(get_curx(profiler)). Profiling is interrupted."
  end

end

function solve_optcache(profiler::ProfilerState)
  optcache = get_solver_state(profiler)
  pars_guess = get_parscache(profiler)
  idx = get_idx(profiler)

  fill_x_reduced!(optcache.reinit_cache.u0, pars_guess, idx)
  set_x_fixed!(optcache.reinit_cache.p, pars_guess[idx])

  return solve!(optcache)
end


function profiler_step!(profiler::ProfilerState, method::IntegrationProfiler)
  
  integrator = get_solver_state(profiler)
  SciMLBase.step!(integrator)

  if SciMLBase.successful_retcode(integrator.sol.retcode)
    idx = get_idx(profiler)

    profiler.curpars .= view(integrator.u, 1:length(integrator.u)-1)
    profiler.curx = integrator.u[idx]
    profiler.curobj = compute_optf(get_optprob(profiler), get_curpars(profiler))
    profiler.solver_retcode = integrator.sol.retcode
    profiler.numiter += 1
  else
    @warn "Solver returned $(integrator.sol.retcode) retcode at profile point x = $(get_curx(profiler)). Profiling is interrupted."
  end

end
