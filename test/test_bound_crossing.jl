using LikelihoodProfiler, OptimizationLBFGSB, ForwardDiff, Test

# The nuisance optimum follows x[2] = 2x[1], leaving a pure quadratic profile.
precision = 1800.0
optf = OptimizationFunction(
  (x, p) -> precision * ((x[1] - 0.7)^2 + (x[2] - 2x[1])^2) / 2,
  AutoForwardDiff()
)
optprob = OptimizationProblem(optf, [0.7, 1.4]; lb=[0.65, 0.0], ub=[0.75, 2.0])
method = OptimizationProfiler(optimizer=LBFGSB(), stepper=FixedStep(initial_step=0.02))

for threshold in (1.92, 3.0)
  @testset "threshold = $threshold" begin
    problem = ProfileLikelihoodProblem(optprob, optprob.u0; idxs=1, threshold)
    curve = only(solve(problem, method; maxiters=100).profiles)
    @test first(curve.x) == 0.65
    @test last(curve.x) == 0.75
    @test curve.x[end-1] ≈ 0.74
    @test all(x -> 0.65 <= x <= 0.75, curve.x)
    @test all(pars -> isapprox(pars[2], 2pars[1]; atol=1e-6), curve.pars)
    interval = endpoints(curve)
    if threshold < precision * 0.05^2 / 2
      @test curve.obj[end-1] < curve.obj_level < curve.obj[end]
      @test retcodes(curve) == (left=:Identifiable, right=:Identifiable)
      @test !isnothing(interval.left) && !isnothing(interval.right)
      if !isnothing(interval.left) && !isnothing(interval.right)
        width = sqrt(2threshold / precision)
        @test interval.left ≈ 0.7 - width atol=5e-4
        @test interval.right ≈ 0.7 + width atol=5e-4
        # Linear interpolation of the known quadratic on the final bracket.
        f_lo, f_hi = precision .* (0.74 - 0.7, 0.75 - 0.7).^2 ./ 2
        expected = 0.74 + 0.01 * (threshold - f_lo) / (f_hi - f_lo)
        @test interval.right ≈ expected atol=1e-6
      end
    else
      @test retcodes(curve) == (left=:NonIdentifiable, right=:NonIdentifiable)
      @test interval == (left=nothing, right=nothing)
    end
  end
end
