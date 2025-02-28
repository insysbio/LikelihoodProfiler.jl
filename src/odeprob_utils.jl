
function solver_init(plprob::PLProblem, method::IntegrationProfiler)
  odeprob = build_scimlprob(plprob, method)

  # If reoptimize is requested, then create an optimization problem, create an
  # optimizer state, and register a callback.

  return SciMLBase.init(odeprob, get_integrator(method); get_integrator_opts(method)...)
end

function solver_reinit!(solver_state::SciMLBase.AbstractODEIntegrator, plprob::PLProblem, method::IntegrationProfiler, idx, dir, profile_bound)
  optpars = get_optpars(plprob)
  x0 = optpars[idx]
  
  # FIXME Bad workaround. We need to update the ODEIntegrator with new tspan and dir.
  odeprob = solver_state.sol.prob
  odeprob2 = remake(odeprob; u0=[optpars; zero(eltype(optpars))], tspan=(x0, profile_bound))
  @show odeprob2.tspan
  @show odeprob2.u0
  #=
  solver_state.tdir = dir
  SciMLBase.reinit!(solver_state, [optpars;0.0]; t0=x0, tf=profile_bound)
  =#

  # update p values
  set_gamma!(odeprob2.p, -get_gamma(method))
  set_idx!(odeprob2.p, idx)
  set_x_fixed!(odeprob2.p, 1.0)

  return SciMLBase.init(odeprob2, get_integrator(method); get_integrator_opts(method)...)
end

#=
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
  set_gamma!(sciml_prob.p, -get_gamma(sciml_prob.p))
  set_idx!(sciml_prob.p, idx)
  set_x_fixed!(sciml_prob.p, 1.0)

  # If reoptimize is requested, then create an optimization problem, create an
  # optimizer state, and register a callback.
  callback = nothing
  if method.reoptimize
    sciml_prob_opt = build_optprob_reduced(plprob.optprob, plprob.optpars)
    solver_state_opt = solver_init(sciml_prob_opt, plprob, method, idx, dir, profile_bound)
    condition(u, t, integrator) = true
    function affect!(integrator)
      set_x_fixed!(solver_state_opt.reinit_cache.p, integrator.u[idx])
      solver_state_opt.reinit_cache.u0 = integrator.u[1:end-1][1:end .!= idx]
      sol = solve!(solver_state_opt)
      for i in 1:length(integrator.u)-1
        i == idx && continue
        integrator.u[i] = sol[i - (i>idx)]
      end
    end
    callback = DiscreteCallback(condition, affect!)
  end
  
  return SciMLBase.init(
    sciml_prob, get_integrator(method);
    get_integrator_opts(method)...,
    callback=callback
  )
end
=#

function build_scimlprob(plprob::PLProblem, method::IntegrationProfiler)
  optprob = get_optprob(plprob)
  optpars = get_optpars(plprob)
  lp = length(optpars)
  optf = OptimizationBase.instantiate_function(optprob.f, optpars, optprob.f.adtype, optprob.p; g=true, h=true)

  odef = build_odefunc(optf, optpars, Val(get_matrix_type(method)))

  gamma = get_gamma(method)
  xspan = (optpars[1], Inf)
  p = FixedParamCache(gamma, 1, 1.0, gamma)

  return ODEProblem(odef, zeros(lp+1), xspan, p)
end

function build_odefunc(optf::OptimizationFunction, optpars, ::Val{:identity})
  lp = length(optpars)
  rhs_vec = similar(optpars)

  function ode_func(dz, z, p, x)
    idx = get_idx(p)
    gamma = get_gamma(p)

    grad! = optf.grad
    grad!(rhs_vec, view(z, 1:lp))
    dz[1:lp] .= .- gamma .* rhs_vec
    dz[idx] = one(dz[idx])
    dz[end] = rhs_vec[idx] + dz[idx]
  end
end

function build_odefunc(optf::OptimizationFunction, optpars, ::Val{:fisher})
  lp = length(optpars)
  T = eltype(optpars)
  lhs_mat = zeros(T, lp, lp)
  rhs_vec = similar(optpars)

  function ode_func(dz, z, p, x)
    #=
    - Fisher information, definition:

        I_ij = 𝔼_Θ (∂ log L / ∂ Θi) (∂ log L / ∂ Θj)

    - We have access to L and ∇L:

        ∂ log L / ∂ Θ  =  (∂ L / ∂ Θ) / L
                              ^
                              ∇L

        I_ij = 𝔼_Θ (∂ L / ∂ Θi) (∂ L / ∂ Θj) / (L^2)

    - Compute I as follows:

        I = 𝔼_[Θ=Θ0] (∇L ∇L.T) / (L^2)

    =#
    # Todo for Sasha: do not use pinv

    idx = get_idx(p)
    gamma = get_gamma(p)

    grad! = optf.grad
    grad!(rhs_vec, view(z, 1:lp))
    # Todo for Sasha: write the correct formula
    rhs_vec = 1 ./ rhs_vec
    lhs_mat = rhs_vec[1:lp] * rhs_vec[1:lp]'

    e_i = zero(z)[1:lp]'
    e_i[idx] = 1
    lhs = [
      lhs_mat   e_i'
      e_i       0
    ]
    rhs = vcat(.- gamma .* (1 ./ rhs_vec), 1)

    dz .= pinv(lhs) * rhs
  end
end

function build_odefunc(optf::OptimizationFunction, optpars, ::Val{:hessian})
  lp = length(optpars)
  T = eltype(optpars)
  lhs_mat = zeros(T, lp, lp)
  rhs_vec = similar(optpars)

 function ode_func(dz, z, p, x)
    #=
    Assume Θ2 is fixed. Solve the following linear system

    ∂^2 L / ∂^2 Θ1      ∂^2 L / ∂ Θ1 ∂ Θ2   0           dΘ1/dC        0
    ∂^2 L / ∂ Θ1 ∂ Θ2   ∂^2 L / ∂^2 Θ2      1     *     dΘ2/dC   =    0
    0                   1                   0           dλ/dC         1

    We note that dΘ2/dC = 1 and eliminate it from the system.
    =#

    idx = get_idx(p)

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
      lhs_mat[i, end] = zero(T)
    end
    lhs_mat[idx, end] = one(T)

    fill_x_full!(dz, pinv(lhs_mat)*rhs_vec, idx, 1.0)
  end
end

