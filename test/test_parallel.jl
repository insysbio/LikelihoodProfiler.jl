
using Test, Base.Threads, Distributed
@everywhere using LikelihoodProfiler, Optimization, OrdinaryDiffEq, CICOBase


@info "Running with $(nthreads()) threads and $(nprocs()) processes"

rosenbrock(x,p) = (1.0 - x[1])^2 + 100.0*(x[2] - x[1]^2)^2

x0 = [1., 1.]
optf = OptimizationFunction(rosenbrock, AutoForwardDiff())
optprob = OptimizationProblem(optf, x0)
#sol = solve(optprob, Optimization.LBFGS())
#sol_u = sol.u

plprob = PLProblem(optprob, sol.u, (-5.,5.); threshold = 1.0)

for method in [
  IntegrationProfiler(integrator = Tsit5()),
  OptimizationProfiler(optimizer = Optimization.LBFGS(), stepper = FixedStep(; initial_step=0.1)),
  CICOProfiler(scan_tol = 1e-4)
]
  
  idxs = vcat([1 for _ in 1:500], [2 for _ in 1:500])
  @time sol_serial = profile(plprob, method, idxs=idxs)
  @time sol_threads = profile(plprob, method, idxs=idxs, parallel_type=:threads)
  
  @test [sol_serial.profiles[i].endpoints for i in 1:length(sol_serial)] == [sol_threads.profiles[i].endpoints for i in 1:length(sol_threads)]
  @test [sol_serial.profiles[i].x for i in 1:length(sol_serial)] == [sol_threads.profiles[i].x for i in 1:length(sol_threads)]

end
