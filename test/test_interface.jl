using LikelihoodProfiler, Test, Optimization

# ParameterProfile tests
optprob1 = OptimizationProblem((x,p)->x[1]^2+x[2]^2, [1., 2.0])

@test_throws DimensionMismatch PLProblem(optprob1, [0.,0, 0.]) 
@test_throws ArgumentError PLProblem(optprob1, [0.,0],  [(-5,4), (-2,1)]; threshold = -1)
@test_throws DimensionMismatch PLProblem(optprob1, [0.,0], [(-5,4), (-2,1), (-3, 20)])

@test_throws ArgumentError profile(PLProblem(optprob1, [0.,0]), OptimizationProfiler(optimizer=Optimization.LBFGS(), stepper = FixedStep(; initial_step=0.01))) 
@test_throws ArgumentError profile(PLProblem(optprob1, [0.,0], [(-Inf,4), (-2,1)]), OptimizationProfiler(optimizer=Optimization.LBFGS(), stepper = FixedStep(; initial_step=0.01)), idxs=[1])

plprob = PLProblem(optprob1, [0.,0], (-5,4))
plprob2 = remake(plprob, profile_range = [nothing, (-5,4)] )
@test_throws ArgumentError profile(plprob2, OptimizationProfiler(optimizer=Optimization.LBFGS(), stepper = FixedStep(; initial_step=0.01)), idxs=[1])

# FunctionProfile tests
#=
optf1 = OptimizationFunction((x,p)->x[1]+x[2])
funprob = PLProblem(optprob1, [0.,0], optf1, (-2,1))

@test_throws ArgumentError PLProblem(optprob1, [0.,0], optf1, (5,4)) 
@test_throws ArgumentError remake(funprob, optpars=zeros(4))
=#
