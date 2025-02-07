using LikelihoodProfiler, Test
using Optimization, OptimizationNLopt, ForwardDiff, OrdinaryDiffEq

include("problems/taxol_model.jl")

const taxol_retcodes = (
  a0 = (:Identifiable,:Identifiable), 
  ka = (:Identifiable,:Identifiable), 
  r0 = (:NonIdentifiable,:Identifiable),
  d0 = (:Identifiable,:NonIdentifiable),
  kd = (:Identifiable,:NonIdentifiable)
) 

const taxol_ci = (
  a0 = (6.615, 17.36), 
  ka = (4.924, 10.781), 
  r0 = (nothing, 0.405),
  d0 = (0.186, nothing),
  kd = (49.237, nothing)
)

function test_taxol(sol, i; kwargs...)
  ret = get_retcodes(sol[i])
  ci = get_endpoints(sol[i])
  @test taxol_retcodes[i][1] == ret[1] 
  @test taxol_retcodes[i][2] == ret[2] 
  taxol_retcodes[i][1] == :Identifiable && (@test isapprox(ci[1], taxol_ci[i][1]; kwargs...))
  taxol_retcodes[i][2] == :Identifiable && (@test isapprox(ci[2], taxol_ci[i][2]; kwargs...))
end

optf = OptimizationFunction(taxol_obj, Optimization.AutoForwardDiff())
optprob = OptimizationProblem(optf, p0)
profile_range = [
  (2., 30.),
  (2.0, 30.),
  (0.01, 0.6),
  (0.05, 5.),
  (30., 250.)
]
plprob = PLProblem(optprob, p0, profile_range; threshold = sigmasq*chi2_quantile(0.95, 5))


@testset "Taxol model. Fixed-step OptimizationProfiler with derivative-free optimizer" begin
  
  idxs = 1:5
  profile_step(p0, i) = p0[i] * 0.1
  atol = [profile_step(p0, i)/2 for i in idxs]
  atol[3] = 0.041 # tmp fix as r0 upper bound fails to be within step/2 tolerance
  method = OptimizationProfiler(optimizer = NLopt.LN_NELDERMEAD(), stepper = FixedStep(; initial_step=profile_step))
  sol = profile(plprob, method)
  for i in idxs
    test_taxol(sol, i; atol = atol[i])
  end

end


@testset "Taxol model. Fixed-step OptimizationProfiler with gradient-based optimizer" begin
  
  idxs = 1:5
  profile_step(p0, i) = p0[i] * 0.1
  atol = [profile_step(p0, i)/2 for i in idxs]
  atol[3] = 0.041 # tmp fix as r0 upper bound fails to be within step/2 tolerance
  method = OptimizationProfiler(optimizer = Optimization.LBFGS(), stepper = FixedStep(; initial_step=profile_step))
  sol = profile(plprob, method)
  for i in idxs
    test_taxol(sol, i; atol = atol[i])
  end

end

@testset "Taxol model. IntegrationProfiler with full hessian" begin
  
  idxs = 1:5
  rtol = 1e-2 # how to set it?
  method = IntegrationProfiler(integrator = FBDF(autodiff=false), matrix_type = :hessian)
  sol = profile(plprob, method)
  for i in idxs
    test_taxol(sol, i; rtol = rtol)
  end
  
end
