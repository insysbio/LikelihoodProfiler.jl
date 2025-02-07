using LikelihoodProfiler, Test
using Optimization, OptimizationNLopt, ForwardDiff, OrdinaryDiffEq

include("problems/jak-stat_model.jl")

const jakstat_retcodes = (
  Epo_degradation_BaF3 = (:Identifiable, :Identifiable), 
  k_exp_hetero =         (:NonIdentifiable, :Identifiable), 
  k_exp_homo =           (:Identifiable, :Identifiable),
  k_imp_hetero =         (:Identifiable, :Identifiable),
  k_imp_homo =           (:Identifiable, :NonIdentifiable),
  k_phos =               (:Identifiable, :Identifiable),
  sd_pSTAT5A_rel =       (:Identifiable, :Identifiable),
  sd_pSTAT5B_rel =       (:Identifiable, :Identifiable),
  sd_rSTAT5A_rel =       (:Identifiable, :Identifiable),
) 

const jakstat_ci = (
  Epo_degradation_BaF3 = (-1.729, -1.402), 
  k_exp_hetero =         (nothing, -3.043), 
  k_exp_homo =           (-2.534, -1.923),
  k_imp_hetero =         (-1.903, -1.665),
  k_imp_homo =           (0.156, nothing),
  k_phos =               (4.112, 4.29),
  sd_pSTAT5A_rel =       (0.42, 0.78),
  sd_pSTAT5B_rel =       (0.674, 0.999),
  sd_rSTAT5A_rel =       (0.357, 0.678),
)

function test_jakstat(sol, i; kwargs...)
  ret = get_retcodes(sol[i])
  ci = get_endpoints(sol[i])
  @test jakstat_retcodes[i][1] == ret[1] 
  @test jakstat_retcodes[i][2] == ret[2] 
  jakstat_retcodes[i][1] == :Identifiable && (@test isapprox(ci[1], jakstat_ci[i][1]; kwargs...))
  jakstat_retcodes[i][2] == :Identifiable && (@test isapprox(ci[2], jakstat_ci[i][2]; kwargs...))
end

optpars = log10.(p_best)
optf = OptimizationFunction(jakstat_obj, Optimization.AutoForwardDiff())
optprob = OptimizationProblem(optf, optpars)
profile_range = [
  (-3.,3.),
  (-6.,4.),
  (-3.,3.),
  (-3.,3),
  (0., 6.),
  (2., 6.),
  (0., 2.),
  (0., 2.),
  (0., 2.)
]

plprob = PLProblem(optprob, optpars, profile_range; threshold = chi2_quantile(0.95, 1))

@testset "JAK2-STAT5 model. Fixed-step OptimizationProfiler with derivative-free optimizer" begin
  
  idxs = 1:9
  profile_step(p0, i) = abs(p0[i]) * 0.05
  atol = [profile_step(optpars, i)/2 for i in idxs]
  method = OptimizationProfiler(optimizer = NLopt.LN_NELDERMEAD(), stepper = FixedStep(; initial_step=profile_step))
  sol = profile(plprob, method)
  for i in idxs
    test_jakstat(sol, i; atol = atol[i])
  end

end


# TODO Fix

@testset "JAK2-STAT5 model. Fixed-step OptimizationProfiler with gradient-based optimizer" begin
  
  idxs = 1:9
  profile_step(p0, i) = abs(p0[i]) * 0.05
  atol = [profile_step(optpars, i)/2 for i in idxs]
  method = OptimizationProfiler(optimizer = Optimization.LBFGS(), stepper = FixedStep(; initial_step=profile_step))
  sol = profile(plprob, method)
  for i in idxs
    test_jakstat(sol, i; atol = atol[i])
  end

end


@testset "JAK2-STAT5 model. IntegrationProfiler with full hessian" begin
  
  idxs = 1:9
  rtol = 3e-2 # how to set it?
  method = IntegrationProfiler(integrator = FBDF(autodiff=false), matrix_type = :hessian)
  sol = profile(plprob, method)
  for i in idxs
    test_jakstat(sol, i; rtol = rtol)
  end
  
end
