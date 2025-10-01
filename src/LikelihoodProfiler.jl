module LikelihoodProfiler

using SciMLBase, PreallocationTools
using SimpleUnPack: @unpack
using Reexport
@reexport import SciMLBase: OptimizationFunction, OptimizationProblem, remake, solve, solve!, init
@reexport using DataFrames
using LinearAlgebra, DataInterpolations
using Distributions: quantile, Chisq
using OptimizationBase
using RecipesBase
using Distributed

abstract type AbstractProfileLikelihoodProblem end
abstract type AbstractProfileTarget end
abstract type AbstractSolverCache end
abstract type AbstractProfilerMethod end
abstract type AbstractProfilerStep{S} end

#struct FunctionProfile <: AbstractProfile end
#struct PredictionProfile <: AbstractProfile end

include("problem_interface.jl")
include("profile_methods.jl")
include("caches.jl")
include("solve.jl")
include("profile_solution.jl")
include("utils.jl")
include("optprob_utils.jl")
include("odeprob_utils.jl")
include("profiler_step.jl")
include("plotting.jl")
include("deprecated.jl")

export PLProblem, ProfileLikelihoodProblem, ProfileLikelihoodSolution
export ParameterTarget, FunctionTarget
export profile #deprecated
export FixedStep # LineSearchStep, InterpolationLineSearch
export chi2_quantile
export OptimizationProfiler, IntegrationProfiler, CICOProfiler
export endpoints, stats, retcodes, obj_level

end #module LikelihoodProfiler
