# push!(LOAD_PATH, "Y:\\")

using LikelihoodProfiler, Test

include("./cases_func.jl")

println("Starting tests for get_interval")
@testset "get_interval" begin include("test_get_interval.jl") end

println("Starting tests for get_endpoint")
@testset "get_endpoint" begin include("test_get_endpoint.jl") end

println("Starting tests for errors")
@testset "errors" begin include("test_errors.jl") end

println("Starting tests for profile")
@testset "profile" begin include("test_profile.jl") end

println("Starting tests for CICO")
@testset "CICO_ONE_PASS" begin include("test_cico_one_pass.jl") end

println("Starting tests for CICO different fitting alg")
@testset "CICO_ONE_PASS_ALG" begin include("test_cico_one_pass_algorithms.jl") end

println("Starting tests for LIN_EXTRAPOL")
@testset "LIN_EXTRAPOL" begin include("test_lin_extrapol.jl") end

println("Starting tests for QUADR_EXTRAPOL")
@testset "QUADR_EXTRAPOL" begin include("test_quadr_extrapol.jl") end
