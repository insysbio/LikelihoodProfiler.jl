using LikelihoodProfiler, Test
using OptimizationLBFGSB, OrdinaryDiffEq, CICOBase

include(joinpath(@__DIR__, "../models/SIR/sir_model.jl"))

const sir_retcodes = ((:Identifiable,:Identifiable), (:Identifiable,:Identifiable), (:Identifiable,:Identifiable)) 
const sir_ci = ((0.376, 0.428), (0.222, 0.277), (1/1.458, 1/1.045))

const rtol = 1e-2

function test_sir(sol, i; kwargs...)
  ret = retcodes(sol[i])
  ci = endpoints(sol[i])
  @test sir_retcodes[i][1] == ret[1] 
  @test sir_retcodes[i][2] == ret[2] 
  sir_retcodes[i][1] == :Identifiable && (@test isapprox(ci[1], sir_ci[i][1]; kwargs...))
  sir_retcodes[i][2] == :Identifiable && (@test isapprox(ci[2], sir_ci[i][2]; kwargs...))
end

optf = OptimizationFunction(sir_obj, AutoForwardDiff())
optprob = OptimizationProblem(optf, p0; lb=[1e-3, 1e-3, 1e-3], ub=[1e3, 1e3, 1.e3])
sol = solve(optprob, LBFGSB())

optpars = sol.u 
plprob = ProfileLikelihoodProblem(optprob, optpars; threshold = chi2_quantile(0.95, 3)/2)

@testset "SIR model. Fixed-step OptimizationProfiler with gradient-based optimizer" begin
  
  method = OptimizationProfiler(optimizer = LBFGSB(), stepper = AdaptiveStep())
  sol = solve(plprob, method)
  for i in eachindex(p0)
    test_sir(sol, i; rtol)
  end

end

@testset "SIR model. Fixed-step OptimizationProfiler with gradient-based optimizer" begin

  method = OptimizationProfiler(optimizer = LBFGSB(), stepper = FixedStep())
  sol = solve(plprob, method)
  for i in eachindex(p0)
    test_sir(sol, i; rtol)
  end

end


@testset "SIR model. IntegrationProfiler with full hessian" begin
  
  method = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (reltol=1e-3, abstol=1e-3), matrix_type = :hessian)
  sol = solve(plprob, method)
  for i in eachindex(p0)
    test_sir(sol, i; rtol)
  end
  
end


@testset "SIR model. CICOProfiler" begin
  
  method = CICOProfiler(optimizer = :LN_NELDERMEAD, scan_tol = 1e-4)
  sol = solve(plprob, method)
  for i in eachindex(p0)
    test_sir(sol, i; rtol)
  end
  
end
