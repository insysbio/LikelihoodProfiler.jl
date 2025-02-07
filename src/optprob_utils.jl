#=
function solver_init(plprob::PLProblem, method::OptimizationProfiler)
  optprob_reduced = build_scimlprob(plprob, method)
  return SciMLBase.init(optprob_reduced, get_optimizer(method), get_optimizer_opts(method)...)
end

function solver_reinit!(solver_state::OptimizationBase.OptimizationCache, plprob::PLProblem, method::OptimizationProfiler, idx, dir, profiler_bound)
  optprob = get_optprob(plprob)
  
  if !isnothing(solver_state.lb)
    lb_full = optprob.lb
    lb_reduced = solver_state.lb
    fill_x_reduced!(lb_reduced, lb_full, idx)
  end

  if !isnothing(solver_state.ub)
    ub_full = optprob.ub
    ub_reduced = solver_state.ub
    fill_x_reduced!(ub_reduced, ub_full, idx)
  end

  u0_full = get_optpars(plprob)
  u0_reduced = solver_state.reinit_cache.u0
  fill_x_reduced!(u0_reduced, u0_full, idx)

  set_idx!(solver_state.reinit_cache.p, idx)
  set_x_fixed!(solver_state.reinit_cache.p, u0_full[idx])
  #p0 = FixedParamCache(solver_state.reinit_cache.p, idx, u0_full[idx])
  #SciMLBase.reinit!(solver_state; p=p0)
  return nothing
end
=#

function solver_init(sciml_prob::SciMLBase.OptimizationProblem, 
  plprob::PLProblem, method::OptimizationProfiler, idx, dir, profile_bound)
  
  optprob = get_optprob(plprob)
  u0_full = get_optpars(plprob)

  # update lb values
  if !isnothing(sciml_prob.lb)
    lb_full = optprob.lb
    lb_reduced = sciml_prob.lb
    fill_x_reduced!(lb_reduced, lb_full, idx)
  end

  # update ub values
  if !isnothing(sciml_prob.ub)
    ub_full = optprob.ub
    ub_reduced = sciml_prob.ub
    fill_x_reduced!(ub_reduced, ub_full, idx)
  end

  # update u0 values
  u0_reduced = sciml_prob.u0
  fill_x_reduced!(u0_reduced, u0_full, idx)


  # update p values
  set_idx!(sciml_prob.p, idx)
  set_x_fixed!(sciml_prob.p, u0_full[idx])
  
  return SciMLBase.init(sciml_prob, get_optimizer(method); get_optimizer_opts(method)...)
end

function build_scimlprob(plprob::PLProblem, method::OptimizationProfiler)
  optprob = get_optprob(plprob)
  optpars = get_optpars(plprob)

  #=
  if length(optpars) == 1 
    return remake(optprob, p=FixedParamCache(optprob.p, 1, optpars[1]))
  else
    =#
    return build_optprob_reduced(optprob, optpars)
  #end
end

function build_optprob_reduced(optprob::OptimizationProblem, optpars)
  #remake(optprob, u0=optpars[idx], lb=optpars[idx], ub=optpars[idx], f=optf_reduced)
  lenpars = length(optpars)
  lenreduced = lenpars - 1

  optf_reduced = build_optf_reduced(optprob, optpars)
  
  lb_full = optprob.lb
  if isnothing(lb_full)
    lb_reduced = nothing
  else
    lb_reduced = zeros(eltype(lb_full), lenreduced)
    #=
    for i in 1:idx-1
      lb_reduced[i] = lb_full[i]
    end
    for i in idx+1:lenpars
      lb_reduced[i-1] = lb_full[i]
    end
    =#
  end

  ub_full = optprob.ub
  if isnothing(ub_full)
    ub_reduced = nothing
  else
    ub_reduced = zeros(eltype(ub_full), lenreduced)
    #=
    for i in 1:idx-1
      ub_reduced[i] = ub_full[i]
    end
    for i in idx+1:lenpars
      ub_reduced[i-1] = ub_full[i]
    end
    =#
  end

  u0_full = optprob.u0
  u0_reduced = zeros(eltype(u0_full), lenreduced)
  #=
  for i in 1:idx-1
    u0_reduced[i] = u0_full[i]
  end
  for i in idx+1:lenpars
    u0_reduced[i-1] = u0_full[i]
  end
  =#
  #OptimizationProblem(optf_reduced, u0_reduced, lb_reduced, ub_reduced)
  remake(optprob, f=optf_reduced, u0=u0_reduced, p=FixedParamCache(optprob.p, 1, u0_full[1]), lb=lb_reduced, ub=ub_reduced)
end

#build_optf_reduced(plprob::PLProblem, method) = nothing
#build_optf_reduced(plprob::PLProblem, method::OptimizationProfiler) =
#  build_optf_reduced(get_optprob(plprob), get_optpars(plprob))

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
      #@show H_full
      #@show x_full
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
