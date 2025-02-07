using LikelihoodProfiler, Test

@testset "Problem interface" begin
  include("test_prob_interface.jl")
end

@testset "Analytic functions" begin
  include("test_analytic_funcs.jl")
end

@testset "SIR model" begin
  include("test_sir_model.jl")
end

@testset "Taxol model" begin
  include("test_taxol_model.jl")
end

@testset "JAK2STAT5 model" begin
  include("test_jak-stat_model.jl")
end