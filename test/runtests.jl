
using LikelihoodProfiler, Base.Test

println("Starting tests for ONE_PASS")
@time @testset "Simple loss functions" begin include("simple_loss_func_test_one_pass.jl") end

println("Starting tests for D2D_PLE")
@time @testset "Simple loss functions" begin include("simple_loss_func_test_d2d_ple.jl") end
