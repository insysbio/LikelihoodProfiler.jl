__precompile__(false)

module LikelihoodProfiler

# include
include("get_endpoint.jl")
include("get_interval.jl")
include("cico_one_pass.jl")

# export
export get_right_endpoint,
    get_endpoint,
    scaling,
    unscaling,
    ProfilePoint,
    EndPoint,
    ParamIntervalInput,
    ParamInterval,
    get_interval
end
