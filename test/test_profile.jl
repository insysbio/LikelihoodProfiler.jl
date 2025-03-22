using LikelihoodProfiler, Test
using Optimization, OptimizationNLopt, Plots, OrdinaryDiffEq, ForwardDiff, CICOBase

######################################### PLProblem ##########################################

optf = OptimizationFunction((x,p) -> 5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2, AutoForwardDiff())
optprob = OptimizationProblem(optf, [0.,0.]; lb=[-3,-1], ub=[10,9])

plprob = PLProblem(optprob, [3.,4.], [(-5,20), (-2,15)]; threshold=4.0)

#################################### OptimizationProfiler ####################################

method1 = OptimizationProfiler(optimizer=NLopt.LN_NELDERMEAD(), stepper = FixedStep(; initial_step=0.1))
sol1 = profile(plprob, method1; verbose=true)


#################################### IntegrationProfiler ####################################

method2 = IntegrationProfiler(integrator=FBDF(), integrator_opts=(dtmax=0.1,), matrix_type=:hessian)
sol2 = profile(plprob, method2; verbose=true)

######################################## CICOProfiler #######################################

method3 = CICOProfiler(optimizer=:LN_NELDERMEAD, scan_tol=1e-3)
sol3 = profile(plprob, method3; verbose=true)


######################################### PLProblem w parameters ##########################################

optf = OptimizationFunction((x,p) -> 5.0 + (p[1] - x[1])^2 + p[2]*(x[2] - x[1]^2)^2, AutoForwardDiff())
optprob = OptimizationProblem(optf, [1.,1.], [1.0, 50.0]; lb=[-100., -100. ], ub=[100., 100.])

plprob = PLProblem(optprob, [1.,1.], (-20,20); threshold=4.0)

#################################### OptimizationProfiler ####################################

method1 = OptimizationProfiler(optimizer=Optimization.LBFGS(), stepper = FixedStep(; initial_step=0.01))
sol1 = profile(plprob, method1; verbose=true)

#################################### IntegrationProfiler ####################################

method2 = IntegrationProfiler(integrator=FBDF(autodiff = AutoFiniteDiff()), integrator_opts=(dtmax=0.1,), matrix_type=:hessian)
sol2 = profile(plprob, method2; verbose=true)

######################################## CICOProfiler #######################################

method3 = CICOProfiler(optimizer=:LN_NELDERMEAD, scan_tol=1e-14)
sol3 = profile(plprob, method3; verbose=true)
