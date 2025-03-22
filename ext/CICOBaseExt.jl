module CICOBaseExt

using LikelihoodProfiler, CICOBase

function LikelihoodProfiler.__profile_dir(plprob::PLProblem, method::CICOProfiler, idx::Int, dir::Int; verbose=false, kwargs...)
  
  # TODO add check_prob_method()
  !LikelihoodProfiler.hasthreshold(plprob) && throw(ArgumentError("`CICOProfiler` doesn't support profiling with infinite `threshold`. Use other profile methods."))
  
  verbose && @info "Computing initial values."
  optprob = LikelihoodProfiler.get_optprob(plprob)
  optpars = LikelihoodProfiler.get_optpars(plprob)
  threshold = LikelihoodProfiler.get_threshold(plprob)
  
  x0 = optpars[idx]
  obj0 = LikelihoodProfiler.compute_optf(optprob, optpars)
  obj_level = obj0 + threshold
  
  profile_range = LikelihoodProfiler.get_profile_range(plprob)
  profile_lb, profile_ub = profile_range isa Tuple ? profile_range : profile_range[idx]
  
  optprob_lb = isnothing(optprob.lb) ? -Inf : optprob.lb
  optprob_ub = isnothing(optprob.ub) ? Inf : optprob.ub
  optprob_range = tuple.(optprob_lb, optprob_ub)
  optprob_bounds = optprob_range isa Tuple ? fill(optprob_range, length(optpars)) : optprob_range

  if dir == 1 
    direction = :right
    profile_bound = profile_ub
  else
    direction = :left
    profile_bound = profile_lb
  end

  verbose && @info "Computing $direction-side profile"

  ep = CICOBase.get_endpoint(optpars, idx, x->optprob.f.f(x,optprob.p), :CICO_ONE_PASS, direction;
    loss_crit = obj_level, theta_bounds=optprob_bounds, scan_bound=profile_bound, local_alg=LikelihoodProfiler.get_optimizer(method), scan_tol=LikelihoodProfiler.get_scan_tol(method), silent=!verbose)

  return cico_to_profile_values(plprob, ep, optpars, x0, obj0, obj_level)
end

function cico_to_profile_values(plprob::PLProblem, ep::CICOBase.EndPoint, optpars, x0, obj0, obj_level)

  ep_retcode = cico_deduce_retcode(ep.status)
  ep_val = ep.value

  if isleft(ep)
    pars = ep_retcode == :Identifiable ? [ep.profilePoints[1].params, optpars] : [optpars]
    x = ep_retcode == :Identifiable ? [ep_val, x0] : [x0]
    obj = ep_retcode == :Identifiable ? [ep.profilePoints[1].loss, obj0] : [obj0]
    retcodes = (ep_retcode, :NonIdentifiable)
    endpoints = (ep_val, nothing)
    stats = ((nf=ep.counter,), nothing)
  else
    pars = ep_retcode == :Identifiable ? [optpars, ep.profilePoints[1].params] : [optpars]
    x = ep_retcode == :Identifiable ? [x0, ep_val] : [x0]
    obj = ep_retcode == :Identifiable ? [obj0, ep.profilePoints[1].loss] : [obj0]
    retcodes = (:NonIdentifiable, ep_retcode)
    endpoints = (nothing, ep_val) 
    stats = (nothing, nf=ep.counter)
  end

  return LikelihoodProfiler.ProfileValues(Val(false), plprob, pars, x, obj, obj_level, retcodes, endpoints, stats)
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

isleft(ep::CICOBase.EndPoint) = ep.direction == :left


end #module