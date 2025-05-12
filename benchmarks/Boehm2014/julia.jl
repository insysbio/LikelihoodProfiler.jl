# 1. Imports
using Revise
using PEtab, LikelihoodProfiler, ForwardDiff, Optimization, OptimizationNLopt, OrdinaryDiffEq, Sundials

###

# 2. Load the model
model_name = "Boehm_JProteomeRes2014"
path_yaml = joinpath(@__DIR__, "../../models/", "$model_name/$model_name.yaml")
petab_model = PEtabModel(path_yaml)

# 3. Define optimization problem
osolver = ODESolver(CVODE_BDF(); abstol=1e-6, reltol=1e-6, abstol_adj = 1e-6, reltol_adj = 1e-6)
petab_problem = PEtabODEProblem(petab_model; gradient_method = :ForwardDiff, hessian_method = :ForwardDiff,
  odesolver = osolver, odesolver_gradient = osolver)
optprob = OptimizationProblem(petab_problem)

# 4. Define profile likelihood problem
optpars = petab_problem.xnominal_transformed
optpars = collect(optpars)
plprob = PLProblem(optprob, optpars)

# 5. Profile
parnames = petab_problem.xnames[1:1]
profile_idxs = collect(1:length(parnames))

# PL methods
## OptimizationProfiler
optmeth = OptimizationProfiler(optimizer = NLopt.LD_LBFGS(), stepper = FixedStep(; initial_step=0.01))
sol1 = profile(plprob, optmeth; idxs=profile_idxs, verbose=true)

# 6. Display results
println("CI, method: $optmeth")
for i in 1:length(parnames)
    println(parnames[i], ", CI is ", get_endpoints(sol1[i]))
end

#=
# Simple example for comparison
rosenbrock(x,p) = (1.0 - x[1])^2 + 100.0*(x[2] - x[1]^2)^2
x0 = zeros(2)
optf = OptimizationFunction(rosenbrock, AutoForwardDiff())
optprob = OptimizationProblem(optf, x0)
sol = solve(optprob, Optimization.LBFGS())
sol_u = sol.u
plprob = PLProblem(optprob, sol_u, (-5.,5.))
method = IntegrationProfiler(
    integrator = Tsit5(), 
    integrator_opts = (dtmax=0.3,), 
    matrix_type = :hessian
)
sol1 = profile(plprob, method; idxs=[1], verbose=true)
=#

