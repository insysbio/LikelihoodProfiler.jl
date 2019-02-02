# specific tests for CICO_ONE_PASS

@testset "f_1p" begin
    res0 = [get_right_endpoint(
        [3.0],
        (x::Vector{Float64}) -> x[i],
        (x::Vector{Float64}) -> f_1p(x) - 9.,
        Val(:CICO_ONE_PASS)
    ) for i in 1:1]
    @test isapprox(res0[1][1], 5., atol=1e-2)
    @test res0[1][3] == :BORDER_FOUND_BY_SCAN_TOL
end
