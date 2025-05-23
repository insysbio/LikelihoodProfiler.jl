using LikelihoodProfiler
using Test, Optimization, OptimizationNLopt, ForwardDiff, OrdinaryDiffEq, CICOBase

@warn "Using $(Base.Threads.nthreads()) threads"

rosenbrock(x,p) = (1.0 - x[1])^2 + 100.0*(x[2] - x[1]^2)^2

x0 = [1., 1.]
optf = OptimizationFunction(rosenbrock, AutoForwardDiff())
optprob = OptimizationProblem(optf, x0)
sol = solve(optprob, Optimization.LBFGS())
sol_u = sol.u

plprob = PLProblem(optprob, sol_u, (-5.,5.); threshold = 1.0)

for method in [
      IntegrationProfiler(integrator = Rosenbrock23(autodiff=false)),
      OptimizationProfiler(optimizer = Optimization.LBFGS(), stepper = FixedStep(; initial_step=0.1))
    ]
  
  idxs = vcat([1 for _ in 1:500], [2 for _ in 1:500])
  @time sol1 = profile(plprob, method, idxs=idxs)
  @time sol2 = profile(plprob, method, idxs=idxs, parallel_type=:threads)
  
  @test [sol1.profiles[i].endpoints for i in 1:length(sol1)] == [sol2.profiles[i].endpoints for i in 1:length(sol2)]
  @test [sol1.profiles[i].x for i in 1:length(sol1)] == [sol2.profiles[i].x for i in 1:length(sol2)]

end
