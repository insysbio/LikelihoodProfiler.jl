using LikelihoodProfiler, Test
using OptimizationNLopt, OptimizationLBFGSB, ForwardDiff, OrdinaryDiffEq, CICOBase

include(joinpath(@__DIR__, "../models/SIR/sir_model.jl"))

const sir_retcodes = ((:Identifiable,:Identifiable), (:Identifiable,:Identifiable), (:Identifiable,:Identifiable)) 
const sir_ci = ((0.376, 0.428), (0.222, 0.277), (1.045e-5, 1.458e-5))

function test_sir(sol, i; kwargs...)
  ret = retcodes(sol[i])
  ci = endpoints(sol[i])
  @test sir_retcodes[i][1] == ret[1] 
  @test sir_retcodes[i][2] == ret[2] 
  sir_retcodes[i][1] == :Identifiable && (@test isapprox(ci[1], sir_ci[i][1]; kwargs...))
  sir_retcodes[i][2] == :Identifiable && (@test isapprox(ci[2], sir_ci[i][2]; kwargs...))
end

optf = OptimizationFunction(sir_obj, AutoForwardDiff())
optprob = OptimizationProblem(optf, p0; lb=[1e-5, 1e-5, 1e-7], ub=[1e5, 1e5, 1e5])
#sol = solve(optprob, NLopt.LN_NELDERMEAD())

optpars = [0.3998583528283355, 0.24676816253516404, 1.2460180516141481e-5]
plprob = ProfileLikelihoodProblem(optprob, optpars; threshold = chi2_quantile(0.95, 3)/2)


@testset "SIR model. Fixed-step OptimizationProfiler with derivative-free optimizer" begin
  
  idxs = 1:3
  profile_step(p0, i) = p0[i] * 0.05
  atol = [profile_step(optpars, i)/2 for i in idxs]
  method = OptimizationProfiler(optimizer = NLopt.LN_NELDERMEAD(), stepper = FixedStep(; initial_step=profile_step))
  sol = solve(plprob, method)
  for i in idxs
    test_sir(sol, i; atol = atol[i])
  end

end

#=TODO Fix

@testset "SIR model. Fixed-step OptimizationProfiler with gradient-based optimizer" begin

  idxs = 1:3
  profile_step(p0, i) = p0[i] * 0.05
  atol = [profile_step(optpars, i)/2 for i in idxs]
  method = OptimizationProfiler(optimizer = LBFGSB(), stepper = FixedStep(; initial_step=profile_step))
  sol = solve(plprob, method)
  for i in idxs
    test_sir(sol, i; atol = atol[i])
  end

end
=# 

@testset "SIR model. IntegrationProfiler with full hessian" begin
  
  idxs = 1:3
  rtol = 3e-3 # how to set it?
  method = IntegrationProfiler(integrator = FBDF(autodiff = AutoFiniteDiff()), matrix_type = :hessian)
  sol = solve(plprob, method)
  for i in idxs
    test_sir(sol, i; rtol = rtol)
  end
  
end


@testset "SIR model. CICOProfiler" begin
  
  idxs = 1:3
  atol = 1e-3
  method = CICOProfiler(optimizer = :LN_NELDERMEAD, scan_tol = 1e-8)
  sol = solve(plprob, method)
  for i in idxs
    test_sir(sol, i; atol = atol)
  end
  
end
