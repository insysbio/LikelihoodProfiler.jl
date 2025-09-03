using Optimization, LikelihoodProfiler, OrdinaryDiffEq, Plots, CICOBase
using PEtab

model_name = "Boehm_JProteomeRes2014"
path_yaml = joinpath(@__DIR__, "../../models/", "$model_name/$model_name.yaml")
petab_model = PEtabModel(path_yaml)
petab_problem = PEtabODEProblem(petab_model)

optprob = OptimizationProblem(petab_problem)

init_pars = get_x(petab_problem)
plprob = ProfileLikelihoodProblem(optprob, init_pars)

alg1 = OptimizationProfiler(; optimizer = Optimization.LBFGS(), stepper = LineSearchStep(; linesearch=InterpolationLineSearch(; min_obj_change=1e-4), initial_step=0.02))
alg2 = OptimizationProfiler(; optimizer = Optimization.LBFGS(), stepper = FixedStep(; initial_step=0.02))

sol1 = solve(plprob, alg1)
sol2 = solve(plprob, alg2)


function plot_sol(sol)
  w, h = 1000, 300

  p11 = plot(sol[1] , dpi=400, xguide="log10_Epo_degradation_BaF3", yguide="likelihood", legend=false)
  p12 = plot(sol[2] , dpi=400, xguide="log10_k_exp_hetero", yguide="likelihood", legend=false, title="OptimizationProfiler")
  p13 = plot(sol[3] , dpi=400, xguide="log10_k_exp_homo", yguide="likelihood", legend=:outerright)

  p1 = plot(p11, p12, p13, layout=(1,3), size=(w,h), dpi=400)

  p21 = plot(sol[4] , dpi=400, xguide="log10_k_imp_hetero", yguide="likelihood", legend=false)
  p22 = plot(sol[5] , dpi=400, xguide="log10_k_imp_homo", yguide="likelihood", legend=false, title="IntegrationProfiler")
  p23 = plot(sol[6] , dpi=400, xguide="log10_k_phos", yguide="likelihood", legend=:outerright)

  p2 = plot(p21, p22, p23, layout=(1,3), size=(w,h), dpi=400)

  p31 = plot(sol[7] , dpi=400, xguide="log10_sd_pSTAT5A_rel", yguide="likelihood", legend=false)
  p32 = plot(sol[8] , dpi=400, xguide="log10_sd_pSTAT5B_rel", yguide="likelihood", legend=false, title="CICOProfiler")
  p33 = plot(sol[9] , dpi=400, xguide="log10_sd_rSTAT5A_rel", yguide="likelihood", legend=:outerright)

  p3 = plot(p31, p32, p33, layout=(1,3), size=(w,h), dpi=400)

  plot(p1, p2, p3, layout=(3,1), size=(w,3*h), dpi=400)
end

plot_sol(sol1)
plot_sol(sol2)