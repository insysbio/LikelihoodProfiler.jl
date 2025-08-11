# 1. Imports
using PEtab, LikelihoodProfiler, ForwardDiff, Optimization, OptimizationNLopt, OrdinaryDiffEq, Sundials
using Plots

# 2. Load the model
model_name = "Boehm_JProteomeRes2014"
path_yaml = joinpath(@__DIR__, "../../models/", "$model_name/$model_name.yaml")
petab_model = PEtabModel(path_yaml)

# 3. Define optimization problem
osolver = ODESolver(CVODE_BDF(); abstol_adj = 1e-3, reltol_adj = 1e-6)
petab_problem = PEtabODEProblem(petab_model; gradient_method = :ForwardDiff, hessian_method = :ForwardDiff,
  odesolver = osolver, odesolver_gradient = osolver)
optprob = OptimizationProblem(petab_problem)

# 4. Define profile likelihood problem
optpars = petab_problem.xnominal_transformed
plprob = PLProblem(optprob, optpars)

# 5. Profile
# All parameters
parnames = petab_problem.xnames
profile_idxs = collect(1:length(parnames))

# PL methods
## OptimizationProfiler
optmeth = OptimizationProfiler(optimizer = NLopt.LD_LBFGS(), stepper = FixedStep(; initial_step=(p0,i)->p0[i]*0.005))
sol1 = solve(plprob, optmeth; idxs=profile_idxs, verbose=true)

## IntegrationProfiler
odemeth = IntegrationProfiler(integrator = FBDF(autodiff = AutoFiniteDiff()), integrator_opts = (dtmax = 0.01, ), matrix_type = :hessian)
sol2 = solve(plprob, odemeth; idxs=profile_idxs, verbose=true)

# 6. Display results
println("CI, method: $optmeth")
for i in 1:length(parnames)
    println(parnames[i], ", CI is ", get_endpoints(sol1[i]))
end

println("CI, method: $odemeth")
for i in 1:length(parnames)
    println(parnames[i], ", CI is ", get_endpoints(sol2[i]))
end
