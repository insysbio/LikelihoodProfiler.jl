
"""
    chi2_quantile(α, df=1)
Computes α quantile for Chi2 distribution with df degrees of freedom
"""
chi2_quantile(α, df=1) = quantile(Chisq(df), α)

evaluate_obj(plprob::ProfileLikelihoodProblem, u, p=plprob.optprob.p) = evaluate_f(plprob.optprob.f, u, p)
evaluate_target_f(plprob::ProfileLikelihoodProblem, idx, u, p=SciMLBase.NullParameters()) = evaluate_target_f(plprob.target, idx, u, p)
evaluate_target_f(target::ParameterTarget, idx, u, p=SciMLBase.NullParameters()) = u[idx]
evaluate_target_f(target::FunctionTarget, idx, u, p=SciMLBase.NullParameters()) = evaluate_f(target.fs[idx], u, p)
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
