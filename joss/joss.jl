using PEtab, Optimization, LikelihoodProfiler, OrdinaryDiffEq, Plots

model_name = "Boehm_JProteomeRes2014"
path_yaml = joinpath(@__DIR__, "../models/", "$model_name/$model_name.yaml")
petab_model = PEtabModel(path_yaml)
petab_problem = PEtabODEProblem(petab_model)


optprob = OptimizationProblem(petab_problem; lb, ub)
plprob = PLProblem(optprob, get_x(petab_problem))

alg1 = OptimizationProfiler(; optimizer = Optimization.LBFGS(), stepper = FixedStep(; initial_step=0.07))
alg2 = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.07,), matrix_type = :hessian)
alg3 = CICOProfiler(optimizer = :LN_NELDERMEAD, scan_tol = 1e-10)


sol1 = profile(plprob, alg1)
sol2 = solve(plprob, alg2)
sol3 = solve(plprob, alg3)

w, h = 1000, 300

p11 = plot(sol1[1] , dpi=400, xguide="log10_Epo_degradation_BaF3", yguide="likelihood", legend=false)
p12 = plot(sol1[2] , dpi=400, xguide="log10_k_exp_hetero", yguide="likelihood", legend=false, title="OptimizationProfiler")
p13 = plot(sol1[3] , dpi=400, xguide="log10_k_exp_homo", yguide="likelihood", legend=:outerright)

p1 = plot(p11, p12, p13, layout=(1,3), size=(w,h), dpi=400)

p21 = plot(sol2[1] , dpi=400, xguide="log10_Epo_degradation_BaF3", yguide="likelihood", legend=false)
p22 = plot(sol2[2] , dpi=400, xguide="log10_k_exp_hetero", yguide="likelihood", legend=false, title="IntegrationProfiler")
p23 = plot(sol2[3] , dpi=400, xguide="log10_k_exp_homo", yguide="likelihood", legend=:outerright)

p2 = plot(p21, p22, p23, layout=(1,3), size=(w,h), dpi=400)

p31 = plot(sol3[1] , dpi=400, xguide="log10_Epo_degradation_BaF3", yguide="likelihood", legend=false)
p32 = plot(sol3[2] , dpi=400, xguide="log10_k_exp_hetero", yguide="likelihood", legend=false, title="CICOProfiler")
p33 = plot(sol3[3] , dpi=400, xguide="log10_k_exp_homo", yguide="likelihood", legend=:outerright)

p3 = plot(p31, p32, p33, layout=(1,3), size=(w,h), dpi=400)

p = plot(p1, p2, p3, layout=(3,1), size=(w,3*h), dpi=400)

savefig(p, joinpath(@__DIR__, "profiles.png"))

df = DataFrame(
    parameter = collect(keys(get_x(petab_problem))),
    OptimizationProfiler = [get_endpoints(sol1[i]) for i in 1:length(sol1)],
    IntegrationProfiler = [get_endpoints(sol2[i]) for i in 1:length(sol2)],
    CICOProfiler = [get_endpoints(sol3[i]) for i in 1:length(sol3)]
)

pretty_table(df; backend = Val(:markdown))
