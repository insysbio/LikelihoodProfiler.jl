
function build_odeprob_reduced(plprob::ProfileLikelihoodProblem, method::IntegrationProfiler, idx, dir, profile_range)
  optprob = plprob.optprob
  optpars = plprob.optpars

  u0 = ComponentVector(theta = optpars, lambda = zero(eltype(optpars)))

  odef = build_odefunc_reduced(optprob, u0, Val(get_matrix_type(method)))

  gamma = get_gamma(method)
  xf = (dir == -1) ? profile_range[1] : profile_range[2]
  xspan = (optpars[idx], xf)
  p = FixedParamCache(optprob.p, idx, 1.0, dir*gamma)

  return ODEProblem(odef, u0, xspan, p)
end

function build_odeprob_full(plprob::ProfileLikelihoodProblem, method::IntegrationProfiler, idx, dir, profile_range)
  optprob = plprob.optprob
  optpars = plprob.optpars

  u0 = ComponentVector(theta = optpars, lambda = zero(eltype(optpars)))

  odef = build_odefunc_full(plprob, u0, idx, Val(get_matrix_type(method)))

  gamma = get_gamma(method)
  x0 = evaluate_target_f(plprob.target, idx, optpars)
  xf = (dir == -1) ? profile_range[1] : profile_range[2]
  xspan = (x0, xf)
  p = FixedParamCache(optprob.p, idx, 1.0, dir*gamma)

  return ODEProblem(odef, u0, xspan, p)
end

function build_odefunc_reduced(optprob::OptimizationProblem, u0, ::Val{:identity})
  rhs_vec = similar(u0.theta)

  optf = OptimizationBase.instantiate_function(optprob.f, u0.theta, optprob.f.adtype, optprob.p; g=true, h=false)

  function ode_func(dz, z, p, x)
    idx = get_idx(p)
    gamma = get_gamma(p)

    grad! = optf.grad
    #grad!(rhs_vec, view(z, 1:lp), p.p)
    grad!(rhs_vec, z.theta, p.p)
    dz.theta .= .- gamma .* rhs_vec
    dz[idx] = one(dz[idx])
    dz[end] = rhs_vec[idx] + dz[idx]
  end
end

#= TODO
function build_odefunc_reduced(optprob::OptimizationProblem, optpars, ::Val{:fisher})
  lp = length(optpars)
  T = eltype(optpars)
  lhs_mat = zeros(T, lp, lp)
  rhs_vec = similar(optpars)

  optf = OptimizationBase.instantiate_function(optprob.f, optpars, optprob.f.adtype, optprob.p; g=true, h=false)


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
    #grad!(rhs_vec, view(z, 1:lp), p.p)
    grad!(rhs_vec, z[1:lp], p.p)
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
=#
function build_odefunc_reduced(optprob::OptimizationProblem, u0, ::Val{:hessian})
  lp = length(u0.theta)
  T = eltype(u0.theta)
  lhs_mat = zeros(T, lp, lp)
  rhs_vec = similar(u0.theta)

  optf = OptimizationBase.instantiate_function(optprob.f, u0.theta, optprob.f.adtype, optprob.p; g=false, h=true)

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
    #hess!(lhs_mat, view(z, 1:lp), p.p)
    #hess!(lhs_mat, view(z, 1:lp))
    hess!(lhs_mat, z.theta, p.p)

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

function build_odefunc_full(plprob::ProfileLikelihoodProblem, u0, idx, ::Val{:hessian})
  optprob = plprob.optprob
  lp = length(u0.theta)
  T = eltype(u0.theta)
  
  lhs_mat = zeros(T, lp+1, lp+1)
  rhs_vec = zeros(T, lp+1)
  rhs_vec[end] = 1.0

  L_hess_mat = zeros(T, lp, lp)
  g_hess_mat = zeros(T, lp, lp)
  g_grad_vec = zeros(T, lp)

  profile_f = plprob.target.fs[idx]
  L_optf = OptimizationBase.instantiate_function(optprob.f, u0.theta, optprob.f.adtype, optprob.p; g=false, h=true)
  g_optf = OptimizationBase.instantiate_function(profile_f, u0.theta, profile_f.adtype, optprob.p; g=true, h=true)

  L_hess! = L_optf.hess
  g_grad! = g_optf.grad
  g_hess! = g_optf.hess

  function ode_func(dz, z, p, x)

    L_hess!(L_hess_mat, z.theta, p.p)
    #@show L_hess_mat
    g_hess!(g_hess_mat, z.theta, p.p)
    g_grad!(g_grad_vec, z.theta, p.p)

    # Top-left block: L_hess_mat .+ z[end] * g_hess_mat
    lhs_mat[1:lp, 1:lp] .= L_hess_mat .+ z.lambda * g_hess_mat

    # Top-right block: g_grad_vec (as a column)
    lhs_mat[1:lp, lp+1] .= g_grad_vec
 
    # Bottom-left block: g_grad_vec' (as a row)
    for i in 1:lp
      lhs_mat[lp+1, i] = g_grad_vec[i]
    end

    # Bottom-right element: 0
    lhs_mat[lp+1, lp+1] = 0

    dz .= pinv(lhs_mat) * rhs_vec
  end
end

function build_odefunc_full(plprob::ProfileLikelihoodProblem, u0, idx, ::Val{:identity})
  optprob = plprob.optprob
  lp = length(u0.theta)
  T = eltype(u0.theta)

  L_grad_vec = zeros(T, lp)
  g_grad_vec = zeros(T, lp)
  
  lhs_mat = zeros(T, lp+1, lp+1)
  rhs_vec = zeros(T, lp+1)
  rhs_vec[end] = 1.0


  profile_f = plprob.target.fs[idx]
  L_optf = OptimizationBase.instantiate_function(optprob.f, u0.theta, optprob.f.adtype, optprob.p; g=true, h=false)
  g_optf = OptimizationBase.instantiate_function(profile_f, u0.theta, profile_f.adtype, optprob.p; g=true, h=false)

  L_grad! = L_optf.grad
  g_grad! = g_optf.grad

  for i in 1:lp
    lhs_mat[i,i] = -1.0
  end

  function ode_func(dz, z, p, x)

    L_grad!(L_grad_vec, z.theta, p.p)
    g_grad!(g_grad_vec, z.theta, p.p)

    gamma = get_gamma(p)
    # Top-right block: g_grad_vec (as a column)
    lhs_mat[1:lp, lp+1] .= g_grad_vec
 
    # Bottom-left block: g_grad_vec' (as a row)
    for i in 1:lp
      lhs_mat[lp+1, i] = g_grad_vec[i]
    end

    # Bottom-right element: 0
    lhs_mat[lp+1, lp+1] = 0

    rhs_vec[1:lp] .= -gamma * L_grad_vec
    dz .= pinv(lhs_mat) * rhs_vec
  end
end

