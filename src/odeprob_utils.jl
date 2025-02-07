#=
function solver_init(plprob::PLProblem, method::IntegrationProfiler)
  odeprob = build_scimlprob(plprob, method)
  return SciMLBase.init(odeprob, get_integrator(method); get_integrator_opts(method)...)
end

function solver_reinit!(solver_state::SciMLBase.AbstractODEIntegrator, plprob::PLProblem, method::IntegrationProfiler, idx, dir, profile_bound)
  optpars = get_optpars(plprob)
  x0 = optpars[idx]
  
  odeprob = solver_state.sol.prob
  odeprob2 = remake(odeprob; u0=[optpars;0.0], tspan=(x0, profile_bound))
  set_idx!(odeprob2.p, idx)
  set_x_fixed!(odeprob2.p, 1.0)
  return SciMLBase.init(odeprob2, get_integrator(method); get_integrator_opts(method)...)
  #=
  set_idx!(solver_state.sol.prob.p, idx)
  set_x_fixed!(solver_state.sol.prob.p, 1.0)
  solver_state.tdir = dir
  SciMLBase.reinit!(solver_state, [optpars;0.0]; t0=x0, tf=profile_bound)
  #@show first_tstop(solver_state)
  return nothing
  =#
end
=#

function solver_init(sciml_prob::SciMLBase.AbstractODEProblem, 
  plprob::PLProblem, method::IntegrationProfiler, idx, dir, profile_bound)

  optpars = get_optpars(plprob)
  x0 = optpars[idx]

  # update u0 values
  for i in eachindex(optpars)
    sciml_prob.u0[i] = optpars[i]
  end
  sciml_prob.u0[end] = 0.0

  # update tspan values
  # Warning: Mutation of ODEProblem detected. SciMLBase v2.0 has made ODEProblem temporarily mutable in order to allow for interfacing with EnzymeRules due to a current limitation in the rule system. This change is only intended to be temporary and ODEProblem will return to being a struct in a later non-breaking release. Do not rely on this behavior, use with caution.
  sciml_prob.tspan = (x0, profile_bound)

  # update p values
  set_idx!(sciml_prob.p, idx)
  set_x_fixed!(sciml_prob.p, 1.0)
  
  return SciMLBase.init(sciml_prob, get_integrator(method); get_integrator_opts(method)...)
end


function build_scimlprob(plprob::PLProblem, method::IntegrationProfiler)
  optprob = get_optprob(plprob)
  optpars = get_optpars(plprob)
  lp = length(optpars)
  optf = OptimizationBase.instantiate_function(optprob.f, optpars, optprob.f.adtype, optprob.p; g=true, h=true)

  odef = build_odefunc(optf, optpars, Val(get_matrix_type(method)))

  gamma = get_gamma(method)
  xspan = (optpars[1], Inf)
  p = FixedParamCache(gamma, 1, 1.0)

  return ODEProblem(odef, zeros(lp+1), xspan, p)
end

function build_odefunc(optf::OptimizationFunction, optpars, ::Val{:hessian})
  lp = length(optpars)
  cache_mat = DiffCache(zeros(lp, lp))
  cache_vec = DiffCache(similar(optpars))

  function ode_func(dz, z, p, x)
    lhs_mat = get_tmp(cache_mat, z)
    rhs_vec = get_tmp(cache_vec, z)
    idx = get_idx(p)
    x_fixed = get_x_fixed(p)

    hess! = optf.hess
    hess!(lhs_mat, view(z, 1:lp))

    # move idx column to the right side
    for i in 1:lp
      rhs_vec[i] = -lhs_mat[i, idx]
    end

    # shift idx:end columns and put lambda column the end
    for i in idx:lp-1
      for j in 1:lp
        lhs_mat[j, i] = lhs_mat[j, i+1]
      end
    end
    for i in 1:lp
      lhs_mat[i, end] = 0.0
    end
    lhs_mat[idx, end] = x_fixed

    fill_x_full!(dz, pinv(lhs_mat)*rhs_vec, idx, x_fixed)
  end
end

