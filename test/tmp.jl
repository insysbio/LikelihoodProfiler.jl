using LikelihoodProfiler, Test
using Optimization, OptimizationNLopt, Plots, OrdinaryDiffEq, ForwardDiff, CICOBase

######################################### PLProblem ##########################################

optf = OptimizationFunction((x,p) -> 5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2, AutoForwardDiff())
optprob = OptimizationProblem(optf, [0.,0.]; lb=[-3,-1], ub=[10,9])

plprob = PLProblem(optprob, [3.,4.], [(-5,20), (-2,15)]; threshold=4.0)

# what happens inside of an ODE system

mat = ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}[ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}}(0.0,0.0) ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}}(0.0,0.0); ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}}(0.0,0.0) ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}}(0.0,0.0)]
v = ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}, Float64, 1}[ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}}(3.0,1.0), ForwardDiff.Dual{ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}}(4.0,0.0)]

optf = OptimizationBase.instantiate_function(optprob.f, [0.,0.], optprob.f.adtype, optprob.p; g=true, h=true)

optf.hess(mat, v) # ERROR