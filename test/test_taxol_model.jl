using LikelihoodProfiler, Test
using OptimizationLBFGSB, OptimizationNLopt, ForwardDiff, OrdinaryDiffEq, CICOBase

include(joinpath(@__DIR__, "../models/Taxol/taxol_model.jl"))

const rtol = 1e-1

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
  ret = retcodes(sol[i])
  ci = endpoints(sol[i])
  @test taxol_retcodes[i][1] == ret[1] 
  @test taxol_retcodes[i][2] == ret[2] 
  taxol_retcodes[i][1] == :Identifiable && (@test isapprox(ci[1], taxol_ci[i][1]; kwargs...))
  taxol_retcodes[i][2] == :Identifiable && (@test isapprox(ci[2], taxol_ci[i][2]; kwargs...))
end

lb = [2.0, 2.0, 0.01, 0.05, 30.0]
ub = [30.0, 30.0, 0.6, 10.0, 210.0]

optf = OptimizationFunction(taxol_obj, AutoForwardDiff())
optprob = OptimizationProblem(optf, p0; lb=lb, ub=ub)

plprob = ProfileLikelihoodProblem(optprob, p0; threshold = sigmasq*chi2_quantile(0.95, 5))

#=
FIXME
@testset "Taxol model. Fixed-step OptimizationProfiler with gradient-based optimizer" begin

  idxs = 1:5
  profile_step(p0, i) = p0[i] * 0.1
  atol = [profile_step(p0, i)/2 for i in idxs]
  method = OptimizationProfiler(; optimizer = LBFGSB(), stepper = FixedStep())
  
  #=
  method = OptimizationProfiler(optimizer = NLopt.LN_NELDERMEAD(), 
  stepper = AdaptiveStep(; initial_step = profile_step,
  controller = ObjectiveStepControl(; target_factor=1.25, step_factor=1.5)))
  =#
  sol = solve(plprob, method)
  for i in idxs
    test_taxol(sol, i; atol = atol[i])
  end

end
=#


@testset "Taxol model. Fixed-step OptimizationProfiler with gradient-based optimizer" begin
  
  method = OptimizationProfiler(optimizer = LBFGSB(), stepper = AdaptiveStep())
  sol = solve(plprob, method)
  for i in eachindex(p0)
    test_taxol(sol, i; rtol = rtol)
  end

end


@testset "Taxol model. IntegrationProfiler with full hessian" begin
  
  method = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (reltol=1e-3, abstol=1e-3), matrix_type = :hessian)
  sol = solve(plprob, method)
  for i in eachindex(p0)
    test_taxol(sol, i; rtol = rtol)
  end
  
end

#=
@testset "Taxol model. CICOProfiler" begin
  
  idxs = 1:5
  atol=1e-3
  method = CICOProfiler(optimizer = :LN_NELDERMEAD, scan_tol = 1e-4)
  sol = solve(plprob, method)
  for i in idxs
    test_taxol(sol, i; atol)
  end
  
end
=#
