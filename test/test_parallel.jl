
using Test, Base.Threads, Distributed

@info "Running with $(nthreads()) threads and $(nprocs()) processes"

addprocs(2)

@everywhere using LikelihoodProfiler, Optimization, ForwardDiff, OrdinaryDiffEq, CICOBase

@everywhere rosenbrock(x,p) = (1.0 - x[1])^2 + 100.0*(x[2] - x[1]^2)^2

x0 = [1., 1.]
optf = OptimizationFunction(rosenbrock, AutoForwardDiff())
optprob = OptimizationProblem(optf, x0)

plprob = ProfileLikelihoodProblem(optprob, x0, (-5.,5.); threshold = 1.0)

for method in [
  IntegrationProfiler(integrator = Tsit5()),
  OptimizationProfiler(optimizer = Optimization.LBFGS(), stepper = FixedStep(; initial_step=0.1)),
  CICOProfiler(scan_tol = 1e-4)
]
  
  idxs = vcat([1 for _ in 1:500], [2 for _ in 1:500])
  @time sol_serial = solve(plprob, method, idxs=idxs)
  @time sol_threads = solve(plprob, method, idxs=idxs, parallel_type=:threads)
  @time sol_distributed = solve(plprob, method, idxs=idxs, parallel_type=:distributed)

  
  @test [sol_serial.profiles[i].endpoints for i in 1:length(sol_serial)] == [sol_threads.profiles[i].endpoints for i in 1:length(sol_threads)]
  @test [sol_serial.profiles[i].x for i in 1:length(sol_serial)] == [sol_threads.profiles[i].x for i in 1:length(sol_threads)]
  @test [sol_serial.profiles[i].endpoints for i in 1:length(sol_serial)] == [sol_distributed.profiles[i].endpoints for i in 1:length(sol_distributed)]
  @test [sol_serial.profiles[i].x for i in 1:length(sol_serial)] == [sol_distributed.profiles[i].x for i in 1:length(sol_distributed)]
end
