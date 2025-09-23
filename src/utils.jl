
"""
    chi2_quantile(α, df=1)
Computes α quantile for Chi2 distribution with df degrees of freedom
"""
chi2_quantile(α, df=1) = quantile(Chisq(df), α)

evaluate_obj(plprob::ProfileLikelihoodProblem, u, p=plprob.optprob.p) = evaluate_f(plprob.optprob.f, u, p)
evaluate_target_f(plprob::ProfileLikelihoodProblem, idx, u, p=SciMLBase.NullParameters()) = evaluate_target_f(plprob.target, idx, u, p)
evaluate_target_f(target::ParameterTarget, idx, u, p=SciMLBase.NullParameters()) = u[idx]
evaluate_targfet_f(target::FunctionTarget, idx, u, p=SciMLBase.NullParameters()) = evaluate_f(target.fs[idx], u, p)
evaluate_f(f::OptimizationFunction, u, p=SciMLBase.NullParameters()) = f(u, p)

evaluate_f(prob::OptimizationProblem, idx, u, p=prob.p) = u[idx]

evaluate_gradf(prob::OptimizationProblem, u, p=prob.p) = evaluate_gradf(prob.f, u, p)
function evaluate_gradf(f::OptimizationFunction, u, p=SciMLBase.NullParameters())
  if !isnothing(f.grad)
    return f.grad(u,p)
  elseif !(f.adtype isa SciMLBase.NoAD)
    OptimizationBase.DifferentiationInterface.gradient(x->f(x,p), f.adtype, u)
  else
    grad_error()
  end
end

evaluate_hessf(prob::OptimizationProblem, u, p=prob.p) = evaluate_hessf(prob.f, u, p)
function evaluate_hessf(f::OptimizationFunction, u, p=SciMLBase.NullParameters())
  if !isnothing(f.hess)
    return f.hess(u,p)
  elseif !(f.adtype isa SciMLBase.NoAD)
    OptimizationBase.DifferentiationInterface.hessian(x->f(x,p), f.adtype, u)
  else
    hess_error()
  end
end 

grad_error() = throw(ArgumentError("Gradient is not provided. Use `OptimizationFunction` either with a valid AD backend https://docs.sciml.ai/Optimization/stable/API/ad/ or a provided 'grad' function."))
hess_error() = throw(ArgumentError("Hessian is not provided. Use `OptimizationFunction` either with a valid AD backend https://docs.sciml.ai/Optimization/stable/API/ad/ or a provided 'hess' function."))


# helper type used in parameter profiling
struct FixedParamCache{P}
  p::P
  idx::Base.RefValue{Int}
  x_fixed::Base.RefValue{Float64}
  gamma::Base.RefValue{Float64}
end

function FixedParamCache(p, idx::Int, x_fixed::Real, gamma::Real)
  # @assert gamma >= 0 # Todo for Sasha: justify
  _p = isnothing(p) ? SciMLBase.NullParameters() : p
  return FixedParamCache{typeof(_p)}(_p, Ref(idx), Ref(float(x_fixed)), Ref(float(gamma)))
end

get_p(p::FixedParamCache) = p.p
get_idx(p::FixedParamCache) = p.idx[]
get_x_fixed(p::FixedParamCache) = p.x_fixed[]
get_gamma(p::FixedParamCache) = p.gamma[]
set_p!(p::FixedParamCache{P}, new_p::P) where {P} = p.p = new_p
set_idx!(p::FixedParamCache, idx::Int) = p.idx[] = idx
set_x_fixed!(p::FixedParamCache, x_fixed::Float64) = p.x_fixed[] = x_fixed
set_gamma!(p::FixedParamCache, gamma::Float64) = p.gamma[] = gamma

fill_x_full!(x_full::AbstractVector{T}, x_reduced::AbstractVector{T}, p::FixedParamCache) where T<:Number = 
  fill_x_full!(x_full, x_reduced, get_idx(p), get_x_fixed(p))

function fill_x_full!(x_full::AbstractVector{T}, x_reduced::AbstractVector{T}, idx::Int, x_fixed) where T<:Number
  for i in 1:idx-1
    x_full[i] = x_reduced[i]
  end
  x_full[idx] = T(x_fixed)
  for i in idx+1:length(x_full)
    x_full[i] = x_reduced[i-1]
  end
end

function fill_x_reduced!(x_reduced::AbstractVector{T}, x_full::AbstractVector{T}, idx::Int) where T
  for i in 1:idx-1
    x_reduced[i] = x_full[i]
  end
  for i in idx+1:length(x_full)
    x_reduced[i-1] = x_full[i]
  end
end

# fill elements of reduced Matrix
function fill_x_reduced!(x_reduced::AbstractMatrix{T}, x_full::AbstractMatrix{T}, idx::Int) where T
  nr, nc = size(x_full)
  for i in 1:idx-1
    for j in 1:idx-1
      x_reduced[i, j] = x_full[i, j]
    end
    for j in idx+1:nc
      x_reduced[i, j-1] = x_full[i, j]
    end
  end

  for i in idx+1:nr
    for j in 1:idx-1
      x_reduced[i-1, j] = x_full[i, j]
    end
    for j in idx+1:nc
      x_reduced[i-1, j-1] = x_full[i, j]
    end
  end
end

get_solver_stats(solver_state::SciMLBase.AbstractODEIntegrator) = solver_state.sol.stats
get_solver_stats(solver_state::SciMLBase.AbstractOptimizationCache) = SciMLBase.OptimizationStats()

