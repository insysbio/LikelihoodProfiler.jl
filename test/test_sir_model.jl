using LikelihoodProfiler, Test
using Optimization, OptimizationNLopt, ForwardDiff, OrdinaryDiffEq

include("problems/sir_model.jl")

const sir_retcodes = ((:Identifiable,:Identifiable), (:Identifiable,:Identifiable), (:Identifiable,:Identifiable)) 
const sir_ci = ((0.376, 0.428), (0.222, 0.277), (1.045e-5, 1.458e-5))

function test_sir(sol, i; kwargs...)
  ret = get_retcodes(sol[i])
  ci = get_endpoints(sol[i])
  @test sir_retcodes[i][1] == ret[1] 
  @test sir_retcodes[i][2] == ret[2] 
  sir_retcodes[i][1] == :Identifiable && (@test isapprox(ci[1], sir_ci[i][1]; kwargs...))
  sir_retcodes[i][2] == :Identifiable && (@test isapprox(ci[2], sir_ci[i][2]; kwargs...))
end

optf = OptimizationFunction(sir_obj, Optimization.AutoForwardDiff())
optprob = OptimizationProblem(optf, p0)
#sol = solve(optprob, NLopt.LN_NELDERMEAD())

optpars = [0.3998583528283355, 0.24676816253516404, 1.2460180516141481e-5]
plprob = PLProblem(optprob, optpars, (-20.,20.); threshold = chi2_quantile(0.95, 3)/2)


@testset "SIR model. Fixed-step OptimizationProfiler with derivative-free optimizer" begin
  
  idxs = 1:3
  profile_step(p0, i) = p0[i] * 0.05
  atol = [profile_step(optpars, i)/2 for i in idxs]
  method = OptimizationProfiler(optimizer = NLopt.LN_NELDERMEAD(), stepper = FixedStep(; initial_step=profile_step))
  sol = profile(plprob, method)
  for i in idxs
    test_sir(sol, i; atol = atol[i])
  end

end


# TODO fix ODE Solver and AutoDiff
# Warning: At t=1.1711631575978938e-15, dt was forced below floating point epsilon 1.9721522630525295e-31, and step error estimate = 34.92360393939161. Aborting. 
#=
@testset "SIR model. Fixed-step OptimizationProfiler with gradient-based optimizer" begin

  idxs = 1:3
  profile_step(p0, i) = p0[i] * 0.05
  atol = [profile_step(optpars, i)/2 for i in idxs]
  method = OptimizationProfiler(optimizer = Optimization.LBFGS(), stepper = FixedStep(; initial_step=profile_step))
  sol = profile(plprob, method)
  for i in idxs
    test_sir(sol, i; atol = atol[i])
  end

end
=#

@testset "SIR model. IntegrationProfiler with full hessian" begin
  
  idxs = 1:3
  rtol = 3e-3 # how to set it?
  method = IntegrationProfiler(integrator = FBDF(autodiff=false), matrix_type = :hessian)
  sol = profile(plprob, method)
  for i in idxs
    test_sir(sol, i; rtol = rtol)
  end
  
end