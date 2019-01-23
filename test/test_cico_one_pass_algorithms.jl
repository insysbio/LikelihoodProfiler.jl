# :LN_COBYLA
@testset "f_3p_1im restricted :LN_COBYLA" begin
    res0 = [get_right_endpoint(
        [3., 8., 2.1],
        i,
        (x::Vector{Float64}) -> f_3p_1im(x) - 9.,
        Val(:CICO_ONE_PASS);
        scan_bound = 4.,
        local_alg = :LN_COBYLA
    ) for i in 1:3]
    @test res0[1][3] == :SCAN_BOUND_REACHED
    @test res0[2][3] == :SCAN_BOUND_REACHED
    @test res0[3][3] == :SCAN_BOUND_REACHED
end

@testset "f_3p_1im_dep :LN_COBYLA" begin
    res0 = [get_right_endpoint(
        [3., 2., 2.1],
        i,
        (x::Vector{Float64}) -> f_3p_1im_dep(x) - 9.,
        Val(:CICO_ONE_PASS);
        local_alg = :LN_COBYLA
    ) for i in 1:3]
    @test isapprox(res0[1][1], 5., atol=1e-2)
    @test res0[1][3] == :BORDER_FOUND_BY_SCAN_TOL
    @test isapprox(res0[2][1], 2.0+2.0*sqrt(2.), atol=1e-2)
    @test res0[2][3] == :BORDER_FOUND_BY_SCAN_TOL
    @test res0[3][3] == :SCAN_BOUND_REACHED
end

@testset "f_3p_1im_dep scan_bound :LN_COBYLA" begin
    res0 = [get_right_endpoint(
        [3., 2., 2.1],
        i,
        (x::Vector{Float64}) -> f_3p_1im_dep(x) - 9.,
        Val(:CICO_ONE_PASS);
        scan_bound = [4.,10.,10.][i],
        local_alg = :LN_COBYLA
    ) for i in 1:3]
    @test res0[1][3] == :SCAN_BOUND_REACHED
    @test isapprox(res0[2][1], 2.0+2.0*sqrt(2.), atol=1e-2)
    @test res0[2][3] == :BORDER_FOUND_BY_SCAN_TOL
    @test res0[3][3] == :SCAN_BOUND_REACHED
end

# :LN_PRAXIS
@testset "f_3p_1im restricted :LN_PRAXIS" begin
    res0 = [get_right_endpoint(
        [3., 8., 2.1],
        i,
        (x::Vector{Float64}) -> f_3p_1im(x) - 9.,
        Val(:CICO_ONE_PASS);
        scan_bound = 4.,
        local_alg = :LN_PRAXIS
    ) for i in 1:3]
    @test res0[1][3] == :SCAN_BOUND_REACHED
    @test res0[2][3] == :SCAN_BOUND_REACHED
    @test res0[3][3] == :SCAN_BOUND_REACHED
end

@testset "f_3p_1im_dep :LN_PRAXIS" begin
    res0 = [get_right_endpoint(
        [3., 2., 2.1],
        i,
        (x::Vector{Float64}) -> f_3p_1im_dep(x) - 9.,
        Val(:CICO_ONE_PASS);
        local_alg = :LN_PRAXIS
    ) for i in 1:3]
    @test isapprox(res0[1][1], 5., atol=1e-2)
    @test res0[1][3] == :BORDER_FOUND_BY_SCAN_TOL
    @test isapprox(res0[2][1], 2.0+2.0*sqrt(2.), atol=1e-2)
    @test res0[2][3] == :BORDER_FOUND_BY_SCAN_TOL
    @test res0[3][3] == :SCAN_BOUND_REACHED
end

@testset "f_3p_1im_dep scan_bound :LN_PRAXIS" begin
    res0 = [get_right_endpoint(
        [3., 2., 2.1],
        i,
        (x::Vector{Float64}) -> f_3p_1im_dep(x) - 9.,
        Val(:CICO_ONE_PASS);
        scan_bound = [4.,10.,10.][i],
        local_alg = :LN_PRAXIS
    ) for i in 1:3]
    @test res0[1][3] == :SCAN_BOUND_REACHED
    @test isapprox(res0[2][1], 2.0+2.0*sqrt(2.), atol=1e-2)
    @test res0[2][3] == :BORDER_FOUND_BY_SCAN_TOL
    @test res0[3][3] == :SCAN_BOUND_REACHED
end
