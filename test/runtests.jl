# using Pkg
# Pkg.add("NLopt")
# Pkg.add("RecipesBase")

using LikelihoodProfiler, Test

include("./cases_func.jl")
methods_list = [:CICO_ONE_PASS, :LIN_EXTRAPOL, :QUADR_EXTRAPOL]

println("Starting tests for get_interval")
@testset "get_interval" begin include("test_get_interval.jl") end

println("Starting tests for get_endpoint")
@testset "get_endpoint CICO_ONE_PASS parameters" begin include("test_get_endpoint_parameters.jl") end
@testset "get_endpoint CICO_ONE_PASS function" begin include("test_get_endpoint_scan_func.jl") end
@testset "get_endpoint errors and warnings" begin include("test_get_endpoint_errors.jl") end

println("Starting tests for all all methods of get_right_endpoint")
@testset "get_right_endpoint methods" begin include("test_methods.jl") end
@testset "get_right_endpoint CICO_ONE_PASS" begin include("test_cico_one_pass.jl") end

println("Starting tests for profile")
@testset "profile" begin include("test_profile.jl") end

println("Starting tests for different fitting alg")
@testset "local_alg" begin include("test_algorithms.jl") end 

println("Starting tests for Plot @recipe")
@testset "PLOT_INTERVAL" begin include("test_plots.jl") end

println("Starting tests for loss error")
@testset "LOSS_ERROR" begin include("test_loss_error.jl") end

println("Starting tests for bands")
@testset "Bands" begin include("test_bands.jl") end

println("Starting tests for get_optimal")
@testset "get_optimal(...)" begin include("test_get_optimal.jl") end

# experimental tests

#@testset "testing derivative-free algorithms" begin include("test_deriv_free_algs.jl") end
#@testset "gradient-based algorithms" begin include("test_grad_algs.jl") end

@testset "get_optimal series" begin include("test_get_optimal_series.jl") end
