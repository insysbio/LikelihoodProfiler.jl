__precompile__(false)
"""
Main module for `LikelihoodProfiler.jl`.

Two functions are exported from this module for public use:

- [`get_endpoint`](@ref). Calculates endpoint of confidence interval.
- [`get_interval`](@ref). Calculates confidence interval.

"""
module LikelihoodProfiler

# include
include("structures.jl")
include("get_endpoint.jl")
include("get_interval.jl")
include("cico_one_pass.jl")
include("method_lin_extrapol.jl")
include("profile.jl")

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
    plot2
end
