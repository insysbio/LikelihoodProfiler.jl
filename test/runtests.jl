# push!(LOAD_PATH, "Y:\\")

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

println("Starting tests for get_interval")
@testset "get_endpoint" begin include("test_param_interval.jl") end

println("Starting tests for CICO different fitting alg")
@testset "CICO_ONE_PASS_ALG" begin include("test_cico_one_pass_algorithms.jl") end

println("Starting tests for errors")
@testset "errors" begin include("test_errors.jl") end
