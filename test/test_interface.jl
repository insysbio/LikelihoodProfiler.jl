using LikelihoodProfiler, Test, Optimization

# ParameterTarget tests
@test_throws DimensionMismatch ParameterTarget(; idxs=1:2, lb=[-10], ub=[10,10])
@test_throws ArgumentError ParameterTarget(; idxs=[1], lb=[7], ub=[4])

optprob1 = OptimizationProblem((x,p)->x[1]^2+x[2]^2, [1., 2.0])
t1 = ParameterTarget(; idxs=1:2, lb=[-5, -2], ub=[4, 1])
t2 = ParameterTarget(; idxs=1:2, lb=[-Inf, -2], ub=[4, 1])

@test_throws DimensionMismatch PLProblem(optprob1, [0.,0, 0.], t1) 
@test_throws ArgumentError PLProblem(optprob1, [0.,0],  t1; threshold = -1)
@test_throws ArgumentError PLProblem(optprob1, [0.,0],  t2)


#plprob = PLProblem(optprob1, [0.,0], t1)

# FunctionProfile tests
#=
optf1 = OptimizationFunction((x,p)->x[1]+x[2])
funprob = PLProblem(optprob1, [0.,0], optf1, (-2,1))

@test_throws ArgumentError PLProblem(optprob1, [0.,0], optf1, (5,4)) 
@test_throws ArgumentError remake(funprob, optpars=zeros(4))
=#
