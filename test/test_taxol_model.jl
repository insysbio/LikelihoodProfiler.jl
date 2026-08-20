using LikelihoodProfiler
using Test
using CICOBase
using OptimizationLBFGSB
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqLowOrderRK: BS3, BS5
using ComponentArrays

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
  ci = endpoints(sol[i], xtransform=exp10)
  @test taxol_retcodes[i][1] == ret[1] 
  @test taxol_retcodes[i][2] == ret[2] 
  taxol_retcodes[i][1] == :Identifiable && (@test isapprox(ci[1], taxol_ci[i][1]; kwargs...))
  taxol_retcodes[i][2] == :Identifiable && (@test isapprox(ci[2], taxol_ci[i][2]; kwargs...))
end

lb = log10.([2.0, 2.0, 0.01, 0.05, 30.0])
ub = log10.([30.0, 30.0, 0.6, 5.0, 170.0])

optf = OptimizationFunction(taxol_obj, AutoForwardDiff())
optprob = OptimizationProblem(optf, p0; lb, ub)

plprob = ProfileLikelihoodProblem(optprob, p0; threshold = sigmasq*chi2_quantile(0.95, 5))

# @btime 34.145 s (802637815 allocations: 20.15 GiB)
@testset "Taxol model. Adaptive-step OptimizationProfiler with gradient-based optimizer" begin
  
  method = OptimizationProfiler(optimizer = LBFGSB(),  stepper = AdaptiveStep())
  sol = solve(plprob, method; reoptimize_init=true)
  for i in eachindex(p0)
    test_taxol(sol, i; rtol)
  end

end

# @btime 15.278 s (192604966 allocations: 34.84 GiB)
@testset "Taxol model. IntegrationProfiler with full hessian" begin
  
  method = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (reltol=1e-4, abstol=1e-6), matrix_type = :hessian)
  sol = solve(plprob, method)
  for i in eachindex(p0)
    test_taxol(sol, i; rtol)
  end
  
end

#= too slow for CI
# @btime 106.840 s (2420937871 allocations: 65.91 GiB)
@testset "Taxol model. IntegrationProfiler with identity and reoptimization" begin

  method = IntegrationProfiler(integrator = BS3(), integrator_opts = (reltol=1e-1, abstol=1e-1), matrix_type = :identity, 
    reoptimize=true, optimizer=LBFGSB())
  sol = solve(plprob, method)
  for i in eachindex(p0)
    test_taxol(sol, i; rtol)
  end
  
end
=#

# @btime 51.595 s (1442105880 allocations: 21.61 GiB)
@testset "Taxol model. CICOProfiler" begin
  
  method = CICOProfiler(optimizer = :LN_NELDERMEAD, scan_tol = 1e-4)
  sol = solve(plprob, method)
  for i in eachindex(p0)
    test_taxol(sol, i; rtol)
  end
  
end

