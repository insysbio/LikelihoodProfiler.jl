using LikelihoodProfiler, OptimizationLBFGSB, OrdinaryDiffEqTsit5, CICOBase, PEtab
using Plots

model_name = "Boehm_JProteomeRes2014"
path_yaml = joinpath(@__DIR__, "../models/", "$model_name/$model_name.yaml")

petab_model = PEtabModel(path_yaml)
petab_problem = PEtabODEProblem(petab_model)

optprob = OptimizationProblem(petab_problem)
optpars = get_x(petab_problem)
param_profile_prob = ProfileLikelihoodProblem(optprob, optpars)

#################################### INTEGRATION-BASED PROFILES ######################################

alg_integ_identity = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.5, reltol=1e-3, abstol=1e-4),
  matrix_type = :identity, gamma=0., reoptimize=true, 
  optimizer = LBFGSB(),optimizer_opts=(maxiters=11000,))

alg_integ_hessian = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.5, reltol=1e-3, abstol=1e-4),
  matrix_type = :hessian, reoptimize=false, 
  optimizer = LBFGSB(),optimizer_opts=(maxiters=11000,))

sol_param_identitiy = solve(param_profile_prob, alg_integ_identity, verbose=true)
sol_param_hessian = solve(param_profile_prob, alg_integ_hessian, verbose=true)

plot(sol_param_identitiy)
#sol_param_fisher = solve(param_profile_prob, alg_integ_fisher, verbose=true)

###################################### CICO PROFILES ######################################

alg_cico = CICOProfiler(optimizer = :LN_NELDERMEAD, scan_tol = 1e-10)
sol_cico = solve(param_profile_prob, alg_cico, verbose=true)

###################################### OPTIMIZATION PROFILES ######################################

alg_optim = OptimizationProfiler(optimizer = LBFGSB(), optimizer_opts=(maxiters=11000,))
sol_param_optim = solve(param_profile_prob, alg_optim, verbose=true)


###################################### FUNCTION PROFILES ######################################

function pSTAT5A_rel_obs(x, p, t)
  sol = get_odesol(x, petab_problem)(t)
  specC17 = 0.107
  return (100 * sol[7] + 200 * sol[6] * specC17) / (sol[7] + sol[1] * specC17 + 2 * sol[6] * specC17)
end

pSTAT5A_rel_optf(t) = OptimizationFunction((x, p) -> pSTAT5A_rel_obs(x, p, t), AutoForwardDiff())

times1 = [60.0, 200.0]
func_profile_prob = ProfileLikelihoodProblem(optprob, optpars, [pSTAT5A_rel_optf(t) for t in times1];
  profile_lower=0.0, profile_upper=120.0)


sol_func_identitiy = solve(func_profile_prob, alg_integ_identity, verbose=true)
sol_func_hessian = solve(func_profile_prob, alg_integ_hessian, verbose=true)
sol_func_optim = solve(func_profile_prob, alg_optim, verbose=true)
sol_func_cico = solve(func_profile_prob, alg_cico, verbose=true)

