module CICOExt

using CICO: get_endpoint, EndPoint

function build_scimlprob(plprob::PLProblem, method::CICOProfiler)
  return nothing
end

function __profile_dir(plprob::PLProblem, method::CICOProfiler, sciml_prob, idx::Int, dir::Int; verbose=true, kwargs...)
  
  verbose && @info "Computing initial values."
  optprob = get_optprob(plprob)
  optpars = get_optpars(plprob)
  threshold = get_threshold(plprob)
  
  obj0 = compute_optf(optprob, optpars)
  obj_level = obj0 + threshold
  
  profile_range = get_profile_range(plprob)
  profile_lb, profile_ub = profile_range isa Tuple ? profile_range : profile_range[idx]
  
  optprob_range = tuple.(optprob.lb, optprob.ub)
  optprob_bounds = optprob_range isa Tuple ? fill(optprob_range, length(optpars)) : optprob_range

  if dir == 1 
    direction = :right
    profile_bound = profile_ub
  else
    direction = :left
    profile_bound = profile_lb
  end

  verbose && @info "Computing $direction-side profile"

  ep = get_endpoint(optpars, idx, x->optprob.f.f(x,SciMLBase.NullParameters()), :CICO_ONE_PASS, direction;
    loss_crit = obj_level, theta_bounds=optprob_bounds, scan_bound=profile_bound, local_alg=get_local_alg(method), scan_tol=get_scan_tol(method), silent=!verbose)

  return cico_to_profile_values(plprob, ep, obj_level)
end

function cico_to_profile_values(plprob::PLProblem, ep::EndPoint, obj_level)

  pars = ep.profilePoints[1].params 
  x = ep.value
  obj = ep.profilePoints[1].loss

  if isleft(ep)
    retcodes = (cico_deduce_retcode(ep.status), :NonIdentifiable)
    endpoints = (x, nothing)
    stats = ((nf=ep.counter,), nothing)
  else
    retcodes = (:NonIdentifiable, cico_deduce_retcode(ep.status))
    endpoints = (nothing ,x) 
    stats = (nothing, nf=ep.counter)
  end
  return ProfileValues(Val(false), plprob, pars, x, obj, obj_level, retcodes, endpoints, stats)
end

function cico_deduce_retcode(retcode::Symbol)
  if retcode == :BORDER_FOUND_BY_SCAN_TOL || retcode == :BORDER_FOUND_BY_LOSS_TOL
    return :Identifiable
  elseif retcode == :SCAN_BOUND_REACHED
    return :NonIdentifiable
  elseif retcode == :MAX_ITER_REACHED
    return :MaxIters
  else
    return :Failed
  end
end

isleft(ep::Endpoint) = ep.direction == :left


end #module