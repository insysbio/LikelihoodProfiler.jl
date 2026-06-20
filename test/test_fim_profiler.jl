using Test, Optimization
using LikelihoodProfiler

const CONF_LEVEL = 0.95
const DF = 1

rosenbrock(x, p) = (1.0 - x[1])^2 + 100.0*(x[2] - x[1]^2)^2

x0 = [1.0, 1.0]
F_true = [802.0  -400.0; -400.0  200.0]
grad_true = [0.0, 0.0]

function rosenbrock_grad!(G, x, p)
  G[1] = -2.0 * (1.0 - x[1]) - 400.0 * x[1] * (x[2] - x[1]^2)
  G[2] = 200.0 * (x[2] - x[1]^2)
  return nothing
end

function rosenbrock_hess!(H, x, p)
  H[1, 1] = 1200.0 * x[1]^2 - 400.0 * x[2] + 2.0
  H[1, 2] = -400.0 * x[1]
  H[2, 1] = -400.0 * x[1]
  H[2, 2] = 200.0
  return nothing
end

optf_iip = OptimizationFunction(rosenbrock; grad=rosenbrock_grad!, hess=rosenbrock_hess!)
optprob_iip = OptimizationProblem(optf_iip, x0)
plprob_iip = ProfileLikelihoodProblem(optprob_iip, x0;
  profile_lower=-10.0, profile_upper=10.0,
  conf_level=CONF_LEVEL, df=DF)

@test LikelihoodProfiler.evaluate_gradf(optprob_iip, x0) == grad_true
@test evaluate_FIM(plprob_iip, x0) == F_true

rosenbrock_grad(x, p) = [
  -2.0 * (1.0 - x[1]) - 400.0 * x[1] * (x[2] - x[1]^2),
  200.0 * (x[2] - x[1]^2)
]
rosenbrock_hess(x, p) = [
  1200.0 * x[1]^2 - 400.0 * x[2] + 2.0  -400.0 * x[1]
  -400.0 * x[1]                              200.0
]

optf_oop = OptimizationFunction(rosenbrock; grad=rosenbrock_grad, hess=rosenbrock_hess)
optprob_oop = OptimizationProblem(optf_oop, x0)
plprob_oop = ProfileLikelihoodProblem(optprob_oop, x0;
  profile_lower=-10.0, profile_upper=10.0,
  conf_level=CONF_LEVEL, df=DF)

F_oop = evaluate_FIM(plprob_oop, x0)
@test LikelihoodProfiler.evaluate_gradf(optprob_oop, x0) == grad_true
@test F_oop == F_true

Finv, _ = LikelihoodProfiler._invert_matrix(F_oop, :cholesky)
Finv_true = [0.5 1.0; 1.0 2.005]
@test isapprox(Finv, Finv_true) 

method = QuadraticApproxProfiler()
sol_oop = solve(plprob_oop, method)
sol_iip = solve(plprob_iip, method)
  
@test size(sol_oop) == (2,)
@test size(sol_iip) == (2,)

for i in 1:length(sol_oop)
  ep_true = (left=x0[i]-sqrt(chi2_quantile(CONF_LEVEL, DF)*Finv_true[i,i]), 
            right=x0[i]+sqrt(chi2_quantile(CONF_LEVEL, DF)*Finv_true[i,i]))
  @test isapprox(endpoints(sol_oop[i]).left, ep_true.left)
  @test isapprox(endpoints(sol_oop[i]).right, ep_true.right)
  @test isapprox(endpoints(sol_iip[i]).left, ep_true.left)
  @test isapprox(endpoints(sol_iip[i]).right, ep_true.right)
  @test retcodes(sol_oop[i]) == (left=:Identifiable, right=:Identifiable)
  @test retcodes(sol_iip[i]) == (left=:Identifiable, right=:Identifiable)
  @test isapprox(sol_oop[i].obj[1], obj_level(sol_oop[1]))
  @test isapprox(sol_iip[i].obj[1], obj_level(sol_iip[1]))
  @test length(sol_oop[i].x) == 101
  @test length(sol_iip[i].x) == 101
end

method2 = QuadraticApproxProfiler(cov_factor=2.0)
sol2 = solve(plprob_oop, method2)
for i in 1:length(sol2)
  ep_true2 = (left=x0[i]-sqrt(2*chi2_quantile(CONF_LEVEL, DF)*Finv_true[i,i]),
              right=x0[i]+sqrt(2*chi2_quantile(CONF_LEVEL, DF)*Finv_true[i,i]))
  @test isapprox(endpoints(sol2[i]).left, ep_true2.left)
  @test isapprox(endpoints(sol2[i]).right, ep_true2.right)
  @test isapprox(sol2[i].obj[1], obj_level(sol2[1]))
end

method3 = QuadraticApproxProfiler(resolution=10)
sol3 = solve(plprob_oop, method3)
for i in 1:length(sol3)
  @test length(sol3[i].x) == 21
end
