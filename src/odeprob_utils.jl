
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

function _ldiv_with_pinv_fallback!(x, A, b, Awork)
  copyto!(Awork, A)
  factorization = lu!(Awork; check=false)

  if LinearAlgebra.issuccess(factorization)
    copyto!(x, b)
    ldiv!(factorization, x)
    all(isfinite, x) && return false
  end

  mul!(x, pinv(A), b)
  return true
end

function build_odefunc_reduced(optprob::OptimizationProblem, u0, ::Val{:identity})
  rhs_vec = similar(u0.theta)

  optf = OptimizationBase.instantiate_function(optprob.f, u0.theta, optprob.f.adtype, optprob.p; g=true, h=false)

  function ode_func(dz, z, p, x)
    idx = get_idx(p)
    gamma = get_gamma(p)

    grad! = optf.grad
    grad!(rhs_vec, z.theta, p.p)
    @inbounds for i in eachindex(rhs_vec)
      dz.theta[i] = -gamma * rhs_vec[i]
    end
    dz.theta[idx] = one(eltype(dz.theta))
    dz.lambda = -gamma * rhs_vec[idx] - dz.theta[idx]
    return nothing
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
  free_indices = zeros(Int, lp-1)
  hess_mat = zeros(T, lp, lp)
  free_hess_mat = zeros(T, lp-1, lp-1)
  factor_mat = similar(free_hess_mat)
  cross_hess_vec = zeros(T, lp-1)
  free_derivative = similar(cross_hess_vec)

  optf = OptimizationBase.instantiate_function(optprob.f, u0.theta, optprob.f.adtype, optprob.p; g=false, h=true)

  function ode_func(dz, z, p, x)
    idx = get_idx(p)

    hess! = optf.hess
    hess!(hess_mat, z.theta, p.p)

    if lp == 1
      dz.theta[1] = one(T)
      dz.lambda = -hess_mat[1, 1]
      return nothing
    end

    free_position = 1
    @inbounds for i in 1:lp
      i == idx && continue
      free_indices[free_position] = i
      free_position += 1
    end

    @inbounds for col in eachindex(free_indices)
      source_col = free_indices[col]
      cross_hess_vec[col] = hess_mat[source_col, idx]
      for row in eachindex(free_indices)
        source_row = free_indices[row]
        free_hess_mat[row, col] = hess_mat[source_row, source_col]
      end
    end

    _ldiv_with_pinv_fallback!(free_derivative, free_hess_mat, cross_hess_vec, factor_mat)
    free_derivative .*= -one(T)

    fill_x_full!(dz.theta, free_derivative, idx, one(T))

    lambda_derivative = -hess_mat[idx, idx]
    @inbounds for i in eachindex(free_indices)
      lambda_derivative -= hess_mat[idx, free_indices[i]] * free_derivative[i]
    end
    dz.lambda = lambda_derivative
    return nothing
  end
end

function build_odefunc_full(plprob::ProfileLikelihoodProblem, u0, idx, ::Val{:hessian})
  optprob = plprob.optprob
  lp = length(u0.theta)
  T = eltype(u0.theta)
  
  lhs_mat = zeros(T, lp+1, lp+1)
  factor_mat = similar(lhs_mat)
  rhs_vec = zeros(T, lp+1)
  rhs_vec[end] = one(T)
  solution_vec = similar(rhs_vec)

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
    g_hess!(g_hess_mat, z.theta, p.p)
    g_grad!(g_grad_vec, z.theta, p.p)

    @inbounds for col in 1:lp
      for row in 1:lp
        lhs_mat[row, col] = muladd(z.lambda, g_hess_mat[row, col], L_hess_mat[row, col])
      end
      lhs_mat[col, lp+1] = g_grad_vec[col]
      lhs_mat[lp+1, col] = g_grad_vec[col]
    end

    lhs_mat[lp+1, lp+1] = zero(T)

    _ldiv_with_pinv_fallback!(solution_vec, lhs_mat, rhs_vec, factor_mat)
    copyto!(dz, solution_vec)
    return nothing
  end
end

function build_odefunc_full(plprob::ProfileLikelihoodProblem, u0, idx, ::Val{:identity})
  optprob = plprob.optprob
  lp = length(u0.theta)
  T = eltype(u0.theta)

  L_grad_vec = zeros(T, lp)
  g_grad_vec = zeros(T, lp)

  profile_f = plprob.target.fs[idx]
  L_optf = OptimizationBase.instantiate_function(optprob.f, u0.theta, optprob.f.adtype, optprob.p; g=true, h=false)
  g_optf = OptimizationBase.instantiate_function(profile_f, u0.theta, profile_f.adtype, optprob.p; g=true, h=false)

  L_grad! = L_optf.grad
  g_grad! = g_optf.grad

  function ode_func(dz, z, p, x)

    L_grad!(L_grad_vec, z.theta, p.p)
    g_grad!(g_grad_vec, z.theta, p.p)

    gamma = get_gamma(p)
    target_grad_norm_sq = dot(g_grad_vec, g_grad_vec)
    iszero(target_grad_norm_sq) && throw(ArgumentError("the profile target gradient must be nonzero"))

    alpha = (one(T) + gamma * dot(g_grad_vec, L_grad_vec)) / target_grad_norm_sq
    @inbounds for i in 1:lp
      dz.theta[i] = muladd(alpha, g_grad_vec[i], -gamma * L_grad_vec[i])
    end
    dz.lambda = -alpha
    return nothing
  end
end

