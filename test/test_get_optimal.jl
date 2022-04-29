# tests for `get_optimal` method
# testing only :LN_PRAXIS method for the case f_1p

@testset "testing get_optimal() for :LN_PRAXIS and f_1p" begin
    res0 = get_optimal(
        [4., 1., 1.],
        f_1p;
        silent = true,
        loss_tol = 1e-3,
        scan_tol = 0.,
        local_alg = :LN_PRAXIS
    )
    @test isapprox(res0.loss, 5., atol=1e-3)
    @test res0.ret == :FTOL_REACHED

    res1 = get_optimal(
        [4., 1., 1.],
        f_1p;
        silent = true,
        loss_tol = 0.,
        scan_tol = 1e-3,
        local_alg = :LN_PRAXIS
    )
    @test isapprox(res1.params[1], 3., atol=1e-3)
    @test res1.ret == :XTOL_REACHED

    res2 = get_optimal(
        [4., 1., 1.],
        f_1p;
        silent = true,
        loss_tol = 0.,
        scan_tol = 0.,
        max_iter = 10,
        local_alg = :LN_PRAXIS
    )
    @test res2.ret == :MAXEVAL_REACHED

    res3 = get_optimal(
        [4., 1., 1.],
        f_1p;
        silent = true,
        loss_tol = 1e-3,
        scan_tol = 0.,
        local_alg = :LN_PRAXIS,
        scale = [:log, :log, :log]
    )
    @test isapprox(res3.loss, 5., atol=1e-3)
    @test res3.ret == :FTOL_REACHED
end
