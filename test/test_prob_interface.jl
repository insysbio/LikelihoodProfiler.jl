using LikelihoodProfiler, Test

# ParameterProfile tests
optprob1 = OptimizationProblem((x,p)->x[1]^2+x[2]^2, [1., 2.0])

@test_throws DimensionMismatch PLProblem(optprob1, [0.])
@test_throws ArgumentError PLProblem(optprob1, [0.,0]) 
@test_throws ArgumentError PLProblem(optprob1, [0.,0, 0.]) 
@test_throws ArgumentError PLProblem(optprob1, [0.,0], (-Inf,4)) 
@test_throws ArgumentError PLProblem(optprob1, [0.,2], (-1,1)) 
@test_throws ArgumentError PLProblem(optprob1, [0.,0],  [(-5,4), (-2,1)]; threshold = -1)
@test_throws DimensionMismatch PLProblem(optprob1, [0.,0], [(-5,4), (-2,1), (-3, 20)])

parprob = PLProblem(optprob1, [0.,0], [(-5,4), (-2,1)])
@test_throws ArgumentError remake(parprob, profile_range=(-Inf,-1))

# FunctionProfile tests
#=
optf1 = OptimizationFunction((x,p)->x[1]+x[2])
funprob = PLProblem(optprob1, [0.,0], optf1, (-2,1))

@test_throws ArgumentError PLProblem(optprob1, [0.,0], optf1, (5,4)) 
@test_throws ArgumentError remake(funprob, optpars=zeros(4))
=#
