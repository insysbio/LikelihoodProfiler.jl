module CICOBaseExt

using LikelihoodProfiler, CICOBase

function LikelihoodProfiler.solve(plprob::ProfileLikelihoodProblem, method::CICOProfiler, idx::Int, dir::Int; verbose=false, kwargs...)
  
  # TODO add check_prob_method()
  !LikelihoodProfiler.hasthreshold(plprob) && throw(ArgumentError("`CICOProfiler` doesn't support profiling with infinite `threshold`. Use other profile methods."))
  
  verbose && @info "Computing initial values."
  optprob = plprob.optprob
  optpars = plprob.optpars
  threshold = plprob.threshold

  obj0 = LikelihoodProfiler.evaluate_obj(plprob, optpars)
  obj_level = obj0 + threshold

  profile_lb = LikelihoodProfiler.get_profile_lb(plprob.target, idx)
  profile_ub = LikelihoodProfiler.get_profile_ub(plprob.target, idx)

  optprob_lb = isnothing(optprob.lb) ? -Inf : float(optprob.lb)
  optprob_ub = isnothing(optprob.ub) ? Inf : float(optprob.ub)
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

  target = plprob.target
  if target isa FunctionTarget
    f = (x) -> target.fs[idx](x, optprob.p)
  elseif target isa ParameterTarget
    f = (x) -> x[idx]
  else
    throw(ArgumentError("Unsupported target type: $(typeof(target))"))
  end

  x0 = LikelihoodProfiler.evaluate_target_f(target, idx, optpars)

  ep = CICOBase.get_endpoint(Vector(optpars), f, x->optprob.f.f(x,optprob.p), :CICO_ONE_PASS, direction;
    loss_crit = obj_level, theta_bounds=Vector(optprob_bounds), scan_bound=profile_bound, local_alg=LikelihoodProfiler.get_optimizer(method), scan_tol=LikelihoodProfiler.get_scan_tol(method), silent=!verbose)

  return cico_to_profile_values(plprob, ep, optpars, idx, dir, x0, obj0, obj_level)
end

function cico_to_profile_values(plprob::ProfileLikelihoodProblem, ep::CICOBase.EndPoint, optpars, idx, dir, x0, obj0, obj_level)

  ep_retcode = cico_deduce_retcode(ep.status)
  ep_val = ep.value

  if isleft(ep)
    pars = ep_retcode == :Identifiable ? [ep.profilePoints[1].params, optpars] : [optpars]
    x = ep_retcode == :Identifiable ? [ep_val, x0] : [x0]
    obj = ep_retcode == :Identifiable ? [ep.profilePoints[1].loss, obj0] : [obj0]
    retcodes = (left = ep_retcode, right = :NotStarted)
    endpoints = (left = ep_val, right = nothing)
    stats = (left = LikelihoodProfiler.SciMLBase.OptimizationStats(;fevals=ep.counter), right = nothing)
  else
    pars = ep_retcode == :Identifiable ? [optpars, ep.profilePoints[1].params] : [optpars]
    x = ep_retcode == :Identifiable ? [x0, ep_val] : [x0]
    obj = ep_retcode == :Identifiable ? [obj0, ep.profilePoints[1].loss] : [obj0]
    retcodes = (left = :NotStarted, right = ep_retcode)
    endpoints = (left = nothing, right = ep_val)
    stats = (left = nothing, right = LikelihoodProfiler.SciMLBase.OptimizationStats(;fevals=ep.counter))
  end

  return LikelihoodProfiler.ProfileCurve{typeof(plprob), typeof(idx), eltype(pars), eltype(x), LikelihoodProfiler.SciMLBase.OptimizationStats}(false, plprob, idx, dir, pars, x, obj, obj_level, retcodes, endpoints, stats)
end

function cico_deduce_retcode(retcode::Symbol)
  if retcode == :BORDER_FOUND_BY_SCAN_TOL || retcode == :BORDER_FOUND_BY_LOSS_TOL
    return :Identifiable
  elseif retcode == :SCAN_BOUND_REACHED
    return :NonIdentifiable
  elseif retcode == :MAX_ITER_REACHED
    return :MaxIters
  else
    return :Failure
  end
end

isleft(ep::CICOBase.EndPoint) = ep.direction == :left


end #module
