using LikelihoodProfiler, OptimizationLBFGSB, OrdinaryDiffEqTsit5, PEtab, CICOBase
using Plots

model_name = "Boehm_JProteomeRes2014"
path_yaml = joinpath(@__DIR__, "../models/", "$model_name/$model_name.yaml")

petab_model = PEtabModel(path_yaml)
petab_problem = PEtabODEProblem(petab_model)

optprob = OptimizationProblem(petab_problem)
optpars = Vector(get_x(petab_problem))
param_profile_prob = ProfileLikelihoodProblem(optprob, optpars; idxs=1:3)

alg_integ = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.5, reltol=1e-3, abstol=1e-4),
  matrix_type = :identity, gamma=0., reoptimize=true, 
  optimizer = LBFGSB(),optimizer_opts=(maxiters=10000,))

###################################### PARAMETER PROFILES ######################################

sol_param = solve(param_profile_prob, alg_integ, verbose=true)

alg_cico = CICOProfiler(optimizer = :LN_NELDERMEAD, scan_tol = 1e-10)
sol_cico = solve(param_profile_prob, alg_cico, verbose=true)

w, h = 1000, 300

p11 = plot(sol_param[1] , dpi=400, xguide="log10_Epo_degradation_BaF3", yguide="objective function", legend=false)
p12 = plot(sol_param[2] , dpi=400, xguide="log10_k_exp_hetero", yguide="objective function", legend=false, title="Parameter Profiles")
p13 = plot(sol_param[3] , dpi=400, xguide="log10_k_exp_homo", yguide="objective function", legend=:outerright)

p1 = plot(p11, p12, p13, layout=(1,3), size=(w,h), margins=5Plots.mm, dpi=400)

savefig(p1, joinpath(@__DIR__, "param_profiles.png"))

p21 = plot(sol_cico[1] , dpi=400, xguide="log10_Epo_degradation_BaF3", yguide="objective function", legend=false)
p22 = plot(sol_cico[2] , dpi=400, xguide="log10_k_exp_hetero", yguide="objective function", legend=false, title="CICOProfiler (Confidence Intervals endpoints)")
p23 = plot(sol_cico[3] , dpi=400, xguide="log10_k_exp_homo", yguide="objective function", legend=:outerright)

p2 = plot(p21, p22, p23, layout=(1,3), size=(w,h), margins=5Plots.mm, dpi=400)

savefig(p2, joinpath(@__DIR__, "cico_profiles.png"))

###################################### FUNCTION PROFILES ######################################

function pSTAT5A_rel_obs(x, p, t)
  sol = get_odesol(x, petab_problem)(t)
  specC17 = 0.107
  return (100 * sol[7] + 200 * sol[6] * specC17) / (sol[7] + sol[1] * specC17 + 2 * sol[6] * specC17)
end

pSTAT5A_rel_optf(t) = OptimizationFunction((x, p) -> pSTAT5A_rel_obs(x, p, t), AutoFiniteDiff())

times1 = [2.5, 60.0, 200.0]
func_profile_prob = ProfileLikelihoodProblem(optprob, optpars, [pSTAT5A_rel_optf(t) for t in times1];
  profile_lower=0.0, profile_upper=120.0)

sol_func = solve(func_profile_prob, alg_integ, verbose=true)

w, h = 1000, 300

p31 = plot(sol_func[1] , dpi=400, xguide="pSTAT5A_rel(t=2.5)", yguide="objective function", legend=false)
p32 = plot(sol_func[2] , dpi=400, xguide="pSTAT5A_rel(t=60)", yguide="objective function", legend=false, title="Function Profiles")
p33 = plot(sol_func[3] , dpi=400, xticks=[25.,30.,35.], xguide="pSTAT5A_rel(t=200)", yguide="objective function", legend=:outerright)

p3 = plot(p31, p32, p33, layout=(1,3), size=(w,h), margins=5Plots.mm, dpi=400)

savefig(p3, joinpath(@__DIR__, "func_profiles.png"))

###################################### PREDICTION PROFILES ######################################

times2 = [2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 160.0, 200.0, 240.0]

pred_profile_prob = ProfileLikelihoodProblem(optprob, optpars, [pSTAT5A_rel_optf(t) for t in times2]; profile_lower=0.0, profile_upper=100.0)

sol_pred = solve(pred_profile_prob, alg_integ, verbose=true)

x = [0; times2]
y = [pSTAT5A_rel_obs(get_x(petab_problem), 0, t) for t in [0; times2]]

# profile endpoints
ep_lower = [0.; [endpoints(sol_pred[i]).left for i in eachindex(times2)]]
ep_upper = [0.; [endpoints(sol_pred[i]).right for i in eachindex(times2)]]

# experimental data
pSTAT5A_rel_data = petab_problem.model_info.petab_measurements.measurement[1:16]

p4 = plot(x, y, ribbon = (y .- ep_lower, ep_upper .- y), xguide="time", lw=2.5, yguide="pSTAT5A_rel", label = "simulation with prediction band", title="Prediction Profile", dpi=400)
scatter!(p4, [0; times2], pSTAT5A_rel_data, label="experimental data")

savefig(p4, joinpath(@__DIR__, "pred_profile.png"))


