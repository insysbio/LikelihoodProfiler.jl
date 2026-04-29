using Test, Optimization
using LikelihoodProfiler

rosenbrock(x, p) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

x0 = [1.0, 1.0]
optf = OptimizationFunction(rosenbrock, AutoForwardDiff())
optprob = OptimizationProblem(optf, x0; lb=[-5.0, -5.0], ub=[5.0, 5.0])
plprob = ProfileLikelihoodProblem(optprob, x0; idxs=[1,2], profile_lower=-5.0, profile_upper=5.0)

method = FIMProfiler()
sol = solve(plprob, method)

@test length(sol) == 2
for i in 1:length(sol)
  ep = endpoints(sol[i])
  @test !isnothing(ep.left)
  @test !isnothing(ep.right)
  @test ep.left <= x0[i] <= ep.right
  @test retcodes(sol[i]).left in (:Identifiable, :NonIdentifiable)
  @test retcodes(sol[i]).right in (:Identifiable, :NonIdentifiable)
end

F = resolve_fim(plprob, FIMProfiler())
@test size(F) == (2, 2)
