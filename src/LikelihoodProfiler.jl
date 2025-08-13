module LikelihoodProfiler

using CommonSolve: CommonSolve, init, solve!, solve
using SciMLBase, PreallocationTools
using Reexport
@reexport import SciMLBase: OptimizationFunction, OptimizationProblem, remake
@reexport using DataFrames
using LineSearch
using LinearAlgebra, DataInterpolations
using Distributions
using OptimizationBase
using RecipesBase
using Distributed

abstract type AbstractProfile end
abstract type AbstractSolverCache end

struct ParameterProfile <: AbstractProfile end
#struct FunctionProfile <: AbstractProfile end
#struct PredictionProfile <: AbstractProfile end

include("problem_interface.jl")
include("profile_solution.jl")
include("profile_methods.jl")
include("profiler_state.jl")
include("profile.jl")
include("utils.jl")
include("optprob_utils.jl")
include("odeprob_utils.jl")
include("endpoints.jl")
include("profiler_step.jl")
include("plotting.jl")
include("deprecated.jl")

export PLProblem, PLSolution
export profile, solve
export FixedStep, LineSearchStep, InterpolationLineSearch
export chi2_quantile
export OptimizationProfiler, IntegrationProfiler, CICOProfiler
export get_endpoints, get_stats, get_retcodes, get_obj_level

end #module LikelihoodProfiler
