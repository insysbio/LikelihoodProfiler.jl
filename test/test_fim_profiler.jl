using Test, Optimization
using LikelihoodProfiler, OptimizationLBFGSB, Plots

const CONF_LEVEL = 0.95
const DF = 1

rosenbrock(x, p) = (1.0 - x[1])^2 + 100.0*(x[2] - x[1]^2)^2

x0 = [1.0, 1.0]
optf = OptimizationFunction(rosenbrock, AutoForwardDiff())
optprob = OptimizationProblem(optf, x0)
plprob = ProfileLikelihoodProblem(optprob, x0; 
  profile_lower=-10.0, profile_upper=10.0, 
  conf_level=CONF_LEVEL, df=DF)

F = evaluate_FIM(plprob, x0)
F_true = [802.0  -400.0; -400.0  200.0]
@test F == F_true

Finv, _ = LikelihoodProfiler._invert_matrix(F, :cholesky)
Finv_true = [0.5 1.0; 1.0 2.005]
@test isapprox(Finv, Finv_true)

evaluate_FIM(plprob, x0)  

method = QuadraticApproxProfiler()
sol = solve(plprob, method)

@test size(sol) == (2,)
for i in 1:length(sol)
  ep_true = (left=x0[i]-sqrt(chi2_quantile(CONF_LEVEL, DF)*Finv_true[i,i]), 
            right=x0[i]+sqrt(chi2_quantile(CONF_LEVEL, DF)*Finv_true[i,i]))
  @test isapprox(endpoints(sol[i]).left, ep_true.left)
  @test isapprox(endpoints(sol[i]).right, ep_true.right)
  @test retcodes(sol[i]) == (left=:Identifiable, right=:Identifiable)
  @test isapprox(sol[i].obj[1], obj_level(sol[1]))
  @test length(sol[i].x) == 101
end

meth_opt = OptimizationProfiler(optimizer = LBFGSB(), stepper = FixedStep(; initial_step=0.15))
sol_opt = solve(plprob, meth_opt)

method2 = QuadraticApproxProfiler(cov_factor=2.0)
sol2 = solve(plprob, method2)
for i in 1:length(sol2)
  ep_true2 = (left=x0[i]-sqrt(2*chi2_quantile(CONF_LEVEL, DF)*Finv_true[i,i]),
              right=x0[i]+sqrt(2*chi2_quantile(CONF_LEVEL, DF)*Finv_true[i,i]))
  @test isapprox(endpoints(sol2[i]).left, ep_true2.left)
  @test isapprox(endpoints(sol2[i]).right, ep_true2.right)
  @test isapprox(sol2[i].obj[1], obj_level(sol2[1]))
end

method3 = QuadraticApproxProfiler(resolution=10)
sol3 = solve(plprob, method3)
for i in 1:length(sol3)
  @test length(sol3[i].x) == 21
end
