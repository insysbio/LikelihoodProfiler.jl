# push!(LOAD_PATH, "Y:\\")
# using Pkg
# Pkg.add("NLopt")
# Pkg.add("RecipesBase")

using LikelihoodProfiler, Test

include("./cases_func.jl")
methods_list = [:CICO_ONE_PASS, :LIN_EXTRAPOL, :QUADR_EXTRAPOL]

println("Starting tests for get_interval")
@testset "get_interval" begin include("test_get_interval.jl") end

println("Starting tests for get_endpoint")
@testset "get_endpoint" begin include("test_get_endpoint.jl") end
@testset "get_endpoint CICO_ONE_PASS" begin include("test_get_endpoint.jl") end
@testset "get_endpoint errors and warnings" begin include("test_get_endpoint_errors.jl") end

println("Starting tests for all all methods of get_right_endpoint")
@testset "get_right_endpoint methods" begin include("test_methods.jl") end
@testset "get_right_endpoint CICO_ONE_PASS" begin include("test_cico_one_pass.jl") end
@testset "get_right_endpoint methods errors and warnings" begin include("test_methods_errors.jl") end

println("Starting tests for profile")
@testset "profile" begin include("test_profile.jl") end

println("Starting tests for different fitting alg")
@testset "local_alg" begin include("test_algorithms.jl") end

println("Starting tests for Plot @recipe")
@testset "PLOT_INTERVAL" begin include("test_plots.jl") end

println("Starting tests for loss error")
@testset "LOSS_ERROR" begin include("test_loss_error.jl") end
