
using Test, Base.Threads, Distributed

@info "Running with $(nthreads()) threads and $(nprocs()) processes"

addprocs(2)

@everywhere using LikelihoodProfiler, OptimizationLBFGSB, ForwardDiff, OrdinaryDiffEq, CICOBase

@everywhere rosenbrock(x,p) = (1.0 - x[1])^2 + 100.0*(x[2] - x[1]^2)^2

x0 = [1., 1.]
optf = OptimizationFunction(rosenbrock, AutoForwardDiff())
optprob = OptimizationProblem(optf, x0)

plprob = ProfileLikelihoodProblem(optprob, x0; profile_lower = -5., profile_upper = 5., threshold = 1.0)

for method in [
  IntegrationProfiler(integrator = Tsit5()),
  OptimizationProfiler(optimizer = LBFGSB(), stepper = FixedStep(; initial_step=0.1)),
  CICOProfiler(scan_tol = 1e-4)
]
  

  @time sol_serial = solve(plprob, method)
  @time sol_threads = solve(plprob, method, parallel_type=:threads)
  @time sol_distributed = solve(plprob, method, parallel_type=:distributed)

  @test [endpoints(sol_serial[i]) for i in 1:length(sol_serial)] == [endpoints(sol_threads[i]) for i in 1:length(sol_threads)]
  @test [sol_serial.profiles[i].x for i in 1:length(sol_serial)] == [sol_threads.profiles[i].x for i in 1:length(sol_threads)]
  @test [sol_serial.profiles[i].x for i in 1:length(sol_serial)] == [sol_distributed.profiles[i].x for i in 1:length(sol_distributed)]
  @test [endpoints(sol_serial[i]) for i in 1:length(sol_serial)] == [endpoints(sol_distributed[i]) for i in 1:length(sol_distributed)]

end
