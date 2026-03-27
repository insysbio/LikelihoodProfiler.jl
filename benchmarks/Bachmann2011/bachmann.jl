# The model is described in the paper "Division of labor by dual feedback regulators controls JAK2/STAT5 signaling over broad ligand range." 10.1038/msb.2011.50

using PEtab, LikelihoodProfiler, ForwardDiff, OptimizationNLopt, OptimizationLBFGSB, OrdinaryDiffEq, Sundials

# Load the model
model_name = "Bachmann_MSB2011"
path_yaml = joinpath(@__DIR__, "../../models/", "$model_name/$model_name.yaml")
petab_model = PEtabModel(path_yaml)

# Optimization problem
osolver = ODESolver(CVODE_BDF(); abstol_adj = 1e-3, reltol_adj = 1e-6)
petab_problem = PEtabODEProblem(petab_model; gradient_method = :ForwardDiff, hessian_method = :ForwardDiff)
#  odesolver = osolver, odesolver_gradient = osolver)
optprob = OptimizationProblem(petab_problem)

# Profile likelihood problem
optpars = petab_problem.xnominal_transformed
plprob = ProfileLikelihoodProblem(optprob, optpars; idxs = 3)

# All parameters
parnames = petab_problem.xnames
# Parameters reported in Identifiability analysis section of the paper
reported_pars = [:CISEqc, :CISEqcOE, :CISInh, :CISRNADelay, 
  :CISRNATurn, :CISTurn, :EpoRActJAK2, :EpoRCISInh, 
  :EpoRCISRemove, :JAK2ActEpo, :JAK2EpoRDeaSHP1, 
  :SHP1ActEpoR, :SHP1Dea, :SHP1ProOE, :SOCS3Eqc,
  :SOCS3EqcOE, :SOCS3Inh, :SOCS3RNADelay, :SOCS3RNATurn,
  :SOCS3Turn, :STAT5ActEpoR, :STAT5ActJAK2, :STAT5Exp,
  :STAT5Imp, :init_EpoRJAK2, :init_SHP1, :init_STAT5]
profile_idxs = [findfirst(isequal(p), parnames) for p in reported_pars]

# PL methods
## OptimizationProfiler
optmeth = OptimizationProfiler(optimizer = NLopt.LD_LBFGS(), stepper = FixedStep(; initial_step=(p0,i)->p0[i]*0.005))
sol = solve(plprob, optmeth; idxs=profile_idxs, verbose=true)

## IntegrationProfiler
alg_integ = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.5, reltol=1e-3, abstol=1e-4),
  matrix_type = :identity, gamma=0., reoptimize=true, 
  optimizer = LBFGSB(),optimizer_opts=(maxiters=11000,))

sol = solve(plprob, alg_integ; verbose=true)
