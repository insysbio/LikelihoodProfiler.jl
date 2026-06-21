using OptimizationLBFGSB, LikelihoodProfiler, OrdinaryDiffEq, Plots, CICOBase
using PEtab

model_name = "Boehm_JProteomeRes2014"
path_yaml = joinpath(@__DIR__, "../../models/", "$model_name/$model_name.yaml")
petab_model = PEtabModel(path_yaml)
petab_problem = PEtabODEProblem(petab_model)

optprob = OptimizationProblem(petab_problem)
init_pars = get_x(petab_problem)

################################### PARAMETER PROFILE LIKELIHOOD PROBLEM ###################################

param_profile_prob = ProfileLikelihoodProblem(optprob, init_pars; idxs=1:3)

alg = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.5, reltol=1e-3, abstol=1e-4), 
  matrix_type = :identity, gamma=0., reoptimize=true, optimizer = LBFGSB(), optimizer_opts=(maxiters=10000,))

sol_param = solve(param_profile_prob, alg, verbose=true)

w, h = 1000, 300

p_param_11 = plot(sol_param[1] , dpi=400, xguide="log10_Epo_degradation_BaF3", yguide="likelihood", legend=false)
p_param_12 = plot(sol_param[2] , dpi=400, xguide="log10_k_exp_hetero", yguide="likelihood", legend=false, title="Parameter Profiles")
p_param_13 = plot(sol_param[3] , dpi=400, xticks=[-2.6,-2.2,-1.8], xguide="log10_k_exp_homo", yguide="likelihood", legend=:outerright)

p_param = plot(p_param_11, p_param_12, p_param_13, layout=(1,3), size=(w,h), margins=5Plots.mm, dpi=400)

################################### FUNCTION PROFILE LIKELIHOOD PROBLEM ###################################

obs = [:pSTAT5A_rel, :pSTAT5B_rel, :rSTAT5A_rel]
times = [2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 160.0, 200.0, 240.0]

function pSTAT5A_rel_obs(x, p, t)
  sol = get_odesol(x, petab_problem)(t)
  specC17 = 0.107
  pSTAT5A_rel = (100 * sol[7] + 200 * sol[6] * specC17) / (sol[7] + sol[1] * specC17 + 2 * sol[6] * specC17)
  return pSTAT5A_rel
end

pSTAT5A_rel_optf(t) = OptimizationFunction((x, p) -> pSTAT5A_rel_obs(x, p, t), AutoForwardDiff())

#=
function pSTAT5B_rel_obs(x, p, t)
  sol = get_odesol(x, petab_problem)(t)
  specC17 = 0.107
    STAT5A = sol[1]
  STAT5B = sol[2]
  pApA = sol[6]
  pApB = sol[7]
   pBpB = sol[8]
  pSTAT5B_rel = -(100 * pApB - 200 * pBpB * (specC17 - 1)) / ((STAT5B * (specC17 - 1) - pApB) + 2 * pBpB * (specC17 - 1))
  return pSTAT5B_rel
end

pSTAT5B_rel_optf(t) = OptimizationFunction((x, p) -> pSTAT5B_rel_obs(x, p, t), AutoFiniteDiff())


function rSTAT5A_rel_obs(x, p, t)
  sol = get_odesol(x, petab_problem)(t)
  STAT5A = sol[1]
  STAT5B = sol[2]
  pApA = sol[6]
  pApB = sol[7]
   pBpB = sol[8]

  specC17 = 0.107
  rSTAT5A_rel = (100 * pApB + 100 * STAT5A * specC17 + 200 * pApA * specC17) / (2 * pApB + STAT5A * specC17 + 2 * pApA * specC17 - STAT5B * (specC17 - 1) - 2 * pBpB * (specC17 - 1))
  return rSTAT5A_rel
end

rSTAT5A_rel_optf(t) = OptimizationFunction((x, p) -> rSTAT5A_rel_obs(x, p, t), AutoFiniteDiff())
=#

func_profile_prob = ProfileLikelihoodProblem(optprob, init_pars, [pSTAT5A_rel_optf(t) for t in [2.5, 60, 200]]; profile_lower=0.0, profile_upper=120.0)

sol_func = solve(func_profile_prob, alg, verbose=true)

w, h = 1000, 300

p_func_11 = plot(sol_func[1] , dpi=400, xguide="pSTAT5A_rel(t=2.5)", yguide="likelihood", legend=false)
p_func_12 = plot(sol_func[2] , dpi=400, xguide="pSTAT5A_rel(t=60)", yguide="likelihood", legend=false, title="Function Profiles")
p_func_13 = plot(sol_func[3] , dpi=400, xticks=[25.,30.,35.], xguide="pSTAT5A_rel(t=200)", yguide="likelihood", legend=:outerright)

p_func = plot(p_func_11, p_func_12, p_func_13, layout=(1,3), size=(w,h), margins=5Plots.mm, dpi=400)

################################### PREDICTION PROFILE LIKELIHOOD PROBLEM ###################################

pred_profile_prob = ProfileLikelihoodProblem(optprob, init_pars, [pSTAT5A_rel_optf(t) for t in times]; profile_lower=0.0, profile_upper=100.0)

sol_pred = solve(pred_profile_prob, alg, verbose=true)

x = [0; times]
y = [pSTAT5A_rel_obs(get_x(petab_problem), 0, t) for t in [0; times]]
ep_lower = [0.; [endpoints(sol_pred[i]).left for i in 1:length(times)]]
ep_upper = [0.; [endpoints(sol_pred[i]).right for i in 1:length(times)]]
pSTAT5A_rel_data = petab_problem.model_info.petab_measurements.measurement[1:16]
p_pred = plot(x, y, ribbon = (y .- ep_lower, ep_upper .- y), xguide="time", lw=2.5, yguide="pSTAT5A_rel", label = "simulation with prediction band", title="Prediction Profile", dpi=400)
scatter!(p_pred, [0; times], pSTAT5A_rel_data, label="experimental data")
