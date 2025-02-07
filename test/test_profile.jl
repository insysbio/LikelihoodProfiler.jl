using LikelihoodProfiler, Test, OptimizationNLopt, Plots, OrdinaryDiffEq, ForwardDiff

######################################### PLProblem ##########################################

optf = OptimizationFunction((x,p) -> 5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2, Optimization.AutoForwardDiff())
optprob = OptimizationProblem(optf, [0.,0.])

plprob = PLProblem(optprob, [3.,4.], [(-5,20), (-2,15)]; threshold=4.0)

#################################### OptimizationProfiler ####################################

method1 = OptimizationProfiler(optimizer=NLopt.LN_NELDERMEAD(), step=FixedStep(; initial_step=0.1))
sol1 = profile(plprob, method1; verbose=true)


#################################### IntegrationProfiler ####################################

method2 = IntegrationProfiler(integrator=FBDF(autodiff=false), integrator_opts=(dtmax=0.1,), matrix_type=:hessian)
sol2 = profile(plprob, method2; verbose=true)

######################################## CICOProfiler #######################################

method3 = CICOProfiler(optimizer=:LN_NELDERMEAD, scan_tol=1e-3)
sol3 = profile(plprob, method3; verbose=true)
