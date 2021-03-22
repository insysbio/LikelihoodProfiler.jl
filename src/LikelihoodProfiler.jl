#__precompile__(false)
"""
Main module for `LikelihoodProfiler.jl`.

Four functions are exported from this module for public use:

- [`get_endpoint`](@ref). Calculates endpoint of confidence interval.
- [`get_interval`](@ref). Calculates confidence interval.
- [`profile`](@ref). Generates the profile function based on `loss_func`
- [`update_profile_points!`](@ref). Updates interval by profile points.

"""
module LikelihoodProfiler

using NLopt, ForwardDiff
using LinearAlgebra
import PlotUtils.adapted_grid

# include
include("structures.jl")
include("get_endpoint.jl")
include("get_endpoint_func.jl")
include("get_interval.jl")
include("cico_one_pass.jl")
include("method_lin_extrapol.jl")
include("method_quadr_extrapol.jl")
include("profile.jl")
include("plot_interval.jl")

# export
export get_right_endpoint,
    get_endpoint,
    scaling,
    unscaling,
    ProfilePoint,
    EndPoint,
    ParamIntervalInput,
    ParamInterval,
    get_interval,
    profile,
    update_profile_points!,
    boxing
end
