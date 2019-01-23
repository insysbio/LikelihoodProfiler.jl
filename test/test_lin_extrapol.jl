method = :LIN_EXTRAPOL

@testset "f_2p_1im" begin
    res0 = [get_right_endpoint(
        [3., 1.],
        i,
        (x::Vector{Float64}) -> f_2p_1im(x) - 9.,
        Val(method)
    ) for i in 1:2]
    @test isapprox(res0[1][1], 5., atol=1e-2)
    @test res0[1][3] == :BORDER_FOUND_BY_SCAN_TOL || res0[1][3] == :BORDER_FOUND_BY_LOSS_TOL
    @test res0[2][3] == :SCAN_BOUND_REACHED
end

@testset "f_2p" begin
    res0 = [get_right_endpoint(
        [3., 4.1],
        i,
        (x::Vector{Float64}) -> f_2p(x) - 9.,
        Val(method);
        scan_tol = 1e-6,
        loss_tol = 0.
    ) for i in 1:2]
    @test isapprox(res0[1][1], 5., atol=1e-6)
    @test res0[1][3] == :BORDER_FOUND_BY_SCAN_TOL || res0[1][3] == :BORDER_FOUND_BY_LOSS_TOL
    @test isapprox(res0[2][1], 6., atol=1e-6)
    @test res0[2][3] == :BORDER_FOUND_BY_SCAN_TOL || res0[2][3] == :BORDER_FOUND_BY_LOSS_TOL
end

@testset "f_3p_1im" begin
    res0 = [get_right_endpoint(
        [3., 8., 2.1],
        i,
        (x::Vector{Float64}) -> f_3p_1im(x) - 9.,
        Val(method)
    ) for i in 1:3]
    @test isapprox(res0[1][1], 5., atol=1e-2)
    @test res0[1][3] == :BORDER_FOUND_BY_SCAN_TOL || res0[1][3] == :BORDER_FOUND_BY_LOSS_TOL
    @test res0[2][3] == :SCAN_BOUND_REACHED
    @test res0[3][3] == :SCAN_BOUND_REACHED
end

@testset "f_3p_1im restricted" begin
    res0 = [get_right_endpoint(
        [3., 8., 2.1],
        i,
        (x::Vector{Float64}) -> f_3p_1im(x) - 9.,
        Val(method);
        scan_bound = 4.
    ) for i in 1:3]
    @test res0[2][3] == :SCAN_BOUND_REACHED
    @test res0[2][3] == :SCAN_BOUND_REACHED
    @test res0[3][3] == :SCAN_BOUND_REACHED
end

@testset "f_3p_1im_dep" begin
    res0 = [get_right_endpoint(
        [3., 2., 2.1],
        i,
        (x::Vector{Float64}) -> f_3p_1im_dep(x) - 9.,
        Val(method)
    ) for i in 1:3]
    @test isapprox(res0[1][1], 5., atol=1e-2)
    @test res0[1][3] == :BORDER_FOUND_BY_SCAN_TOL || res0[1][3] == :BORDER_FOUND_BY_LOSS_TOL
    @test isapprox(res0[2][1], 2.0+2.0*sqrt(2.), atol=1e-2)
    @test res0[2][3] == :BORDER_FOUND_BY_SCAN_TOL || res0[2][3] == :BORDER_FOUND_BY_LOSS_TOL
    @test res0[3][3] == :SCAN_BOUND_REACHED
end

@testset "f_3p_1im_dep scan_bound" begin
    res0 = [get_right_endpoint(
        [3., 2., 2.1],
        i,
        (x::Vector{Float64}) -> f_3p_1im_dep(x) - 9.,
        Val(method);
        scan_bound = [4.,10.,10.][i]
    ) for i in 1:3]
    @test res0[1][3] == :SCAN_BOUND_REACHED
    @test isapprox(res0[2][1], 2.0+2.0*sqrt(2.), atol=1e-2)
    @test res0[2][3] == :BORDER_FOUND_BY_SCAN_TOL || res0[2][3] == :BORDER_FOUND_BY_LOSS_TOL
    @test res0[3][3] == :SCAN_BOUND_REACHED
end

@testset "f_3p_1im_dep scan_tol" begin
    res0 = [get_right_endpoint(
        [3., 2., 2.1],
        1,
        (x::Vector{Float64}) -> f_3p_1im_dep(x) - 9.,
        Val(method);
        scan_tol = [1e-2, 1e-4, 1e-6][i]
    ) for i in 1:3]
    @test isapprox(res0[1][1], 5., atol=1e-1)
    @test isapprox(res0[2][1], 5., atol=1e-3)
    @test isapprox(res0[3][1], 5., atol=1e-5)
end

@testset "f_4p_2im" begin
    res0 = [get_right_endpoint(
        [3., 4, 1.1, 10.],
        i,
        (x::Vector{Float64}) -> f_4p_2im(x) - 9.,
        Val(method)
    ) for i in 1:4]
    @test isapprox(res0[1][1], 5., atol=1e-2)
    @test res0[1][3] == :BORDER_FOUND_BY_SCAN_TOL || res0[1][3] == :BORDER_FOUND_BY_LOSS_TOL
    @test isapprox(res0[2][1], 6., atol=1e-2)
    @test res0[2][3] == :BORDER_FOUND_BY_SCAN_TOL || res0[2][3] == :BORDER_FOUND_BY_LOSS_TOL
    @test res0[3][3] == :SCAN_BOUND_REACHED
    @test res0[4][3] == :SCAN_BOUND_REACHED
end

@testset "f_4p_3im" begin
    res0 = [get_right_endpoint(
        [3., 4, 1.1, 10.],
        i,
        (x::Vector{Float64}) -> f_4p_3im(x) - 9.,
        Val(method)
    ) for i in 1:4]

    @test isapprox(res0[1][1], 5., atol=1e-2)
    @test res0[1][3] == :BORDER_FOUND_BY_SCAN_TOL || res0[1][3] == :BORDER_FOUND_BY_LOSS_TOL
    @test res0[2][3] == :SCAN_BOUND_REACHED
    @test res0[3][3] == :SCAN_BOUND_REACHED
    @test res0[4][3] == :SCAN_BOUND_REACHED
end

@testset "f_1p_ex" begin
    res0 = [get_right_endpoint(
        [1.5],
        i,
        (x::Vector{Float64}) -> f_1p_ex(x) - 9.,
        Val(method)
    ) for i in 1:1]

    @test isapprox(res0[1][1], 2. + 1e-8, atol=1e-2)
    @test res0[1][3] == :BORDER_FOUND_BY_SCAN_TOL || res0[1][3] == :BORDER_FOUND_BY_LOSS_TOL
end

@testset "f_5p_3im" begin
    res0 = [get_right_endpoint(
        [3., 0.1, 4, 1.1, 8.],
        i,
        (x::Vector{Float64}) -> f_5p_3im(x) - 9.,
        Val(method)
    ) for i in 1:5]

    @test isapprox(res0[1][1], 5., atol=1e-2)
    @test res0[1][3] == :BORDER_FOUND_BY_SCAN_TOL || res0[1][3] == :BORDER_FOUND_BY_LOSS_TOL
    @test isapprox(res0[2][1], log(3.), atol=1e-2)
    @test res0[2][3] == :BORDER_FOUND_BY_SCAN_TOL || res0[2][3] == :BORDER_FOUND_BY_LOSS_TOL
    @test res0[3][3] == :SCAN_BOUND_REACHED
    @test res0[4][3] == :SCAN_BOUND_REACHED
    @test res0[5][3] == :SCAN_BOUND_REACHED
end

@testset "f_3p_im" begin
    res0 = [get_right_endpoint(
        [3., 0.1, 8.],
        i,
        (x::Vector{Float64}) -> f_3p_im(x) - 9.,
        Val(method)
    ) for i in 1:3]

    @test isapprox(res0[1][1], 5., atol=1e-2)
    @test res0[1][3] == :BORDER_FOUND_BY_SCAN_TOL || res0[1][3] == :BORDER_FOUND_BY_LOSS_TOL
    @test isapprox(res0[2][1], log(3.), atol=1e-2)
    @test res0[2][3] == :BORDER_FOUND_BY_SCAN_TOL || res0[2][3] == :BORDER_FOUND_BY_LOSS_TOL
    @test res0[3][3] == :SCAN_BOUND_REACHED
end
