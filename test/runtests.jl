# push!(LOAD_PATH, "Y:\\")
# ] add NLopt
# using NLopt

using LikelihoodProfiler
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end
include("./cases_func.jl")

println("Starting tests for CICO")
@testset "CICO_ONE_PASS" begin include("test_cico_one_pass.jl") end

println("Starting tests for get_endpoint")
@testset "get_endpoint" begin include("test_get_endpoint.jl") end
