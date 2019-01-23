@testset "wrong theta_init" begin
    @test_throws ArgumentError get_endpoint(
        [-10., 2., 2.1],
        1,
        f_3p_1im_dep,
        :CICO_ONE_PASS;
        loss_crit = 9.
        )
end

@testset "wrong scan_bound" begin
    @test_throws ArgumentError get_endpoint(
        [3., 2., 2.1],
        1,
        f_3p_1im_dep,
        :CICO_ONE_PASS;
        loss_crit = 9.,
        scan_bound = 6.,
        theta_bounds = [(-Inf,5.),(-Inf,Inf),(-Inf,Inf)]
        )
end

@testset "wrong scan_bound and theta_init" begin
    @test_throws ArgumentError get_endpoint(
        [3., 2., 2.1],
        1,
        f_3p_1im_dep,
        :CICO_ONE_PASS;
        loss_crit = 9.,
        scan_bound = 2.
        )
end

@testset "wrong theta_init in log scale" begin
    @test_throws ArgumentError get_endpoint(
        [3., 2., -1],
        1,
        f_3p_1im_dep,
        :CICO_ONE_PASS;
        loss_crit = 9.,
        scale = [:direct,:direct,:log]
        )
end

@testset "wrong theta_bounds in log scale" begin
    @test_throws ArgumentError get_endpoint(
        [3., 2., 2.],
        3,
        f_3p_1im_dep,
        :CICO_ONE_PASS;
        loss_crit = 9.,
        scale = [:direct,:direct,:log],
        theta_bounds = [(-Inf, Inf), (-Inf, Inf), (-5., Inf)]
        )
end
