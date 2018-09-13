
using ParametersIdentification, Test

println("Starting tests")
# @time @testset "Wrong input tests" begin include("wrong_input_test.jl") end
@time @testset "Simple loss functions" begin include("simple_loss_func_test.jl") end
