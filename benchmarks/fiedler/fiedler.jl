# The model is described in the paper "Division of labor by dual feedback regulators controls JAK2/STAT5 signaling over broad ligand range." 10.1038/msb.2011.50

using PEtab, LikelihoodProfiler, ForwardDiff, Optimization, OptimizationNLopt, OrdinaryDiffEq, Sundials

# Load the model
path_yaml = joinpath(@__DIR__, "Fiedler_BMCSystBiol2016", "Fiedler_BMCSystBiol2016.yaml")
petab_model = PEtabModel(path_yaml)

# Optimization problem
osolver = ODESolver(CVODE_BDF(); abstol_adj = 1e-3, reltol_adj = 1e-6)
petab_problem = PEtabODEProblem(petab_model; gradient_method = :ForwardDiff, hessian_method = :ForwardDiff,
  odesolver = osolver, odesolver_gradient = osolver)
optprob = OptimizationProblem(petab_problem)

# Profile likelihood problem
optpars = petab_problem.xnominal_transformed
plprob = PLProblem(optprob, optpars)

# All parameters
parnames = petab_problem.xnames
profile_idxs = collect(1:length(parnames))

# PL methods
## OptimizationProfiler
optmeth = OptimizationProfiler(optimizer = NLopt.LD_LBFGS(), stepper = FixedStep(; initial_step=(p0,i)->p0[i]*0.005))
# sol = profile(plprob, optmeth; idxs=profile_idxs, verbose=true)

## IntegrationProfiler
odemeth = IntegrationProfiler(integrator = FBDF(autodiff = AutoFiniteDiff()), integrator_opts = (dtmax = 0.01, ), matrix_type = :hessian)
sol = profile(plprob, odemeth; idxs=profile_idxs, verbose=true)
