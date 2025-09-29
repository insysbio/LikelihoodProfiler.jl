

function build_optprob_reduced(optprob::OptimizationProblem, optpars, idx)
  lenpars = length(optpars)
  lenreduced = lenpars - 1

  optf_reduced = build_optf_reduced(optprob, optpars)
  
  lb_full = optprob.lb
  if isnothing(lb_full)
    lb_reduced = nothing
  else
    lb_reduced = zeros(eltype(lb_full), lenreduced)
    fill_x_reduced!(lb_reduced, lb_full, idx)
  end

  ub_full = optprob.ub
  if isnothing(ub_full)
    ub_reduced = nothing
  else
    ub_reduced = zeros(eltype(ub_full), lenreduced)
    fill_x_reduced!(ub_reduced, ub_full, idx)
  end
  
  u0_full = optprob.u0
  u0_reduced = zeros(eltype(u0_full), lenreduced)
  fill_x_reduced!(u0_reduced, u0_full, idx)

  remake(optprob, f=optf_reduced, u0=u0_reduced, p=FixedParamCache(optprob.p, idx, u0_full[idx], 1.0), lb=lb_reduced, ub=ub_reduced)
end

function build_optf_reduced(optprob::OptimizationProblem, optpars)
  lenpars = length(optpars)
  xcache = DiffCache(similar(optpars))
  optf_full = optprob.f

  # reduced obj function
  f_full = optf_full.f
  f_reduced = (x_reduced, p) -> begin
    idx = get_idx(p)
    x_fixed = get_x_fixed(p)
    x_full = get_tmp(xcache, x_reduced)
    fill_x_full!(x_full, x_reduced, idx, x_fixed)
    
    return f_full(x_full,p.p)
  end

  # reduced grad function
  # in general we should check if grad, hess, etc have inplace form: e.g. (G, u, p) or G(u, p)
  grad_full! = optf_full.grad
  if isnothing(grad_full!)
    grad_reduced! = nothing
  else
    grad_cache = DiffCache(similar(optpars))
    grad_reduced! = (G_reduced, x_reduced, p) -> begin
      idx = get_idx(p)
      x_fixed = get_x_fixed(p)
      x_full = get_tmp(xcache, x_reduced)
      fill_x_full!(x_full, x_reduced, idx, x_fixed)
      
      G_full = get_tmp(grad_cache, x_full)
      grad_full!(G_full, x_full, p.p)
      fill_x_reduced!(G_reduced, G_full, idx)

      return nothing
    end
  end

  # reduced hess function
  hess_full! = optf_full.hess
  if isnothing(hess_full!)
    hess_reduced! = nothing
  else
    hess_cache = DiffCache(zeros(eltype(optpars), lenpars, lenpars))
    hess_reduced! = (H_reduced, x_reduced, p) -> begin
      idx = get_idx(p)
      x_fixed = get_x_fixed(p)
      x_full = get_tmp(xcache, x_reduced)
      fill_x_full!(x_full, x_reduced, idx, x_fixed)
      
      H_full = get_tmp(hess_cache, x_full)
      hess_full!(H_full, x_full, p.p)
      fill_x_reduced!(H_reduced, H_full, idx)

      return nothing
    end
  end

  #TODO add hessprod, hess_prototype, etc
  
  kwargs = Dict{Symbol,Any}()
  for p in propertynames(optf_full)
    if p ∉ (:f, :adtype, :grad, :hess)
      kwargs[p] = getfield(optf_full, p)
    end
  end
  return OptimizationFunction(f_reduced, optf_full.adtype; grad=grad_reduced!, hess=hess_reduced!, kwargs...)
end


function build_optprob_constrained(plprob::ProfileLikelihoodProblem, idx)
  optprob = plprob.optprob
  optpars = plprob.optpars

  x0 = evaluate_target_f(plprob.target, idx, optpars)
  optf_constrained = build_optf_constrained(plprob, idx)
  return OptimizationProblem(optf_constrained, copy(optpars);
     lb=optprob.lb, ub=optprob.ub, lcons=[x0], ucons=[x0])
end

function build_optf_constrained(plprob::ProfileLikelihoodProblem, idx)
  optprob = plprob.optprob
  cons_optf = get_profile_fs(plprob.target)[idx]

  function cons_f(res, x, p)
    res[1] = cons_optf.f(x, p)
    return nothing
  end
  return OptimizationFunction(optprob.f.f, cons_optf.adtype; grad=optprob.f.grad, hess=optprob.f.hess, cons = cons_f)
end