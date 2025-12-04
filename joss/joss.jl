using OptimizationLBFGSB, LikelihoodProfiler, OrdinaryDiffEqTsit5, Plots, PEtab, CICOBase

model_name = "Boehm_JProteomeRes2014"
path_yaml = joinpath(@__DIR__, "../models/", "$model_name/$model_name.yaml")

petab_model = PEtabModel(path_yaml)
petab_problem = PEtabODEProblem(petab_model)

optprob = OptimizationProblem(petab_problem)
param_profile_prob = ProfileLikelihoodProblem(optprob, get_x(petab_problem); idxs=1:3)

alg = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.5, reltol=1e-3, abstol=1e-4),
  matrix_type = :identity, gamma=0., reoptimize=true, 
  optimizer = LBFGSB(),optimizer_opts=(maxiters=10000,))

sol_param = solve(param_profile_prob, alg, verbose=true)

alg_cico = CICOProfiler(optimizer = :LN_NELDERMEAD, scan_tol = 1e-10)
sol_cico = solve(param_profile_prob, alg_cico, verbose=true)

w, h = 1000, 300

p11 = plot(sol_param[1] , dpi=400, xguide="log10_Epo_degradation_BaF3", yguide="likelihood", legend=false)
p12 = plot(sol_param[2] , dpi=400, xguide="log10_k_exp_hetero", yguide="likelihood", legend=false, title="Parameter Profiles")
p13 = plot(sol_param[3] , dpi=400, xguide="log10_k_exp_homo", yguide="likelihood", legend=:outerright)

p1 = plot(p11, p12, p13, layout=(1,3), size=(w,h), dpi=400)

savefig(p1, joinpath(@__DIR__, "param_profiles.png"))

p21 = plot(sol_cico[1] , dpi=400, xguide="log10_Epo_degradation_BaF3", yguide="likelihood", legend=false)
p22 = plot(sol_cico[2] , dpi=400, xguide="log10_k_exp_hetero", yguide="likelihood", legend=false, title="CICOProfiler (Confidence Intervals endpoints)")
p23 = plot(sol_cico[3] , dpi=400, xguide="log10_k_exp_homo", yguide="likelihood", legend=:outerright)

p2 = plot(p21, p22, p23, layout=(1,3), size=(w,h), dpi=400)

savefig(p2, joinpath(@__DIR__, "cico_profiles.png"))

function pSTAT5A_rel_obs(x, p, t)
  sol = get_odesol(x, petab_problem)(t)
  specC17 = 0.107
  return (100 * sol[7] + 200 * sol[6] * specC17) / (sol[7] + sol[1] * specC17 + 2 * sol[6] * specC17)
end
pSTAT5A_rel_optf(t) = OptimizationFunction((x, p) -> pSTAT5A_rel_obs(x, p, t), AutoFiniteDiff())

times = [2.5, 60.0, 200.0]
func_profile_prob = ProfileLikelihoodProblem(optprob, get_x(petab_problem), [pSTAT5A_rel_optf(t) for t in times];
  profile_lower=0.0, profile_upper=120.0)

sol_func = solve(func_profile_prob, alg)

