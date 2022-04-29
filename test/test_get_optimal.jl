# tests for `get_optimal` method
# testing only the case f_2p
# using LikelihoodProfiler

@testset "testing get_optimal() for :LN_PRAXIS and f_2p, identifiable" begin
    res0 = get_optimal(
        [4., 1.],
        f_2p;
        silent = true,
        loss_tol = 1e-3,
        scan_tol = 0.,
        local_alg = :LN_PRAXIS
    )
    @test isapprox(res0.loss, 5., atol=1e-3)
    @test res0.ret == :FTOL_REACHED

    res1 = get_optimal(
        [4., 1.],
        f_2p;
        silent = true,
        loss_tol = 0.,
        scan_tol = 1e-3,
        local_alg = :LN_PRAXIS
    )
    @test isapprox(res1.params[1], 3., atol=1e-3)
    @test isapprox(res1.params[2], 4., atol=1e-3)
    @test res1.ret == :XTOL_REACHED

    res2 = get_optimal(
        [4., 1.],
        f_2p;
        silent = true,
        loss_tol = 0.,
        scan_tol = 0.,
        max_iter = 10,
        local_alg = :LN_PRAXIS
    )
    @test res2.ret == :MAXEVAL_REACHED

    res3 = get_optimal(
        [4., 1.],
        f_2p;
        silent = true,
        loss_tol = 1e-3,
        scan_tol = 0.,
        local_alg = :LN_PRAXIS,
        scale = [:log, :log]
    )
    @test isapprox(res3.loss, 5., atol=1e-3)
    @test res3.ret == :FTOL_REACHED
end

@testset "testing get_optimal() for :LN_PRAXIS and f_2p, nonidentifiable" begin
    res0 = get_optimal(
        [4., 1., 1.],
        f_2p;
        silent = true,
        loss_tol = 1e-3,
        scan_tol = 0.,
        local_alg = :LN_PRAXIS
    )
    @test isapprox(res0.loss, 5., atol=1e-3)
    @test res0.ret == :FTOL_REACHED

    res1 = get_optimal(
        [4., 1., 1.],
        f_2p;
        silent = true,
        loss_tol = 0.,
        scan_tol = 1e-3,
        local_alg = :LN_PRAXIS
    )
    @test isapprox(res1.params[1], 3., atol=1e-3)
    @test isapprox(res1.params[2], 4., atol=1e-3)
    @test res1.ret == :XTOL_REACHED

    res2 = get_optimal(
        [4., 1., 1.],
        f_2p;
        silent = true,
        loss_tol = 0.,
        scan_tol = 0.,
        max_iter = 10,
        local_alg = :LN_PRAXIS
    )
    @test res2.ret == :MAXEVAL_REACHED

    res3 = get_optimal(
        [4., 1., 1.],
        f_2p;
        silent = true,
        loss_tol = 1e-3,
        scan_tol = 0.,
        local_alg = :LN_PRAXIS,
        scale = [:log, :log, :log]
    )
    @test isapprox(res3.loss, 5., atol=1e-3)
    @test res3.ret == :FTOL_REACHED
end

f_2p(x) = 5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2
grad_2p(x) = [2*(x[1]-3.), 2*(x[2]-4.), 0.]

@testset "testing get_optimal() for :LD_MMA and f_2p, nonidentifiable" begin
    res0 = get_optimal(
        [4., 1., 1.],
        f_2p;
        silent = true,
        loss_tol = 1e-3,
        scan_tol = 0.,
        local_alg = :LD_MMA,
        loss_grad = grad_2p
    )
    @test isapprox(res0.loss, 5., atol=1e-3)
    @test res0.ret == :FTOL_REACHED

    res1 = get_optimal(
        [4., 1., 1.],
        f_2p;
        silent = true,
        loss_tol = 0.,
        scan_tol = 1e-3,
        local_alg = :LD_MMA,
        loss_grad = grad_2p
    )
    @test isapprox(res1.params[1], 3., atol=1e-3)
    @test isapprox(res1.params[2], 4., atol=1e-3)
    @test res1.ret == :XTOL_REACHED

    res2 = get_optimal(
        [4., 1., 1.],
        f_2p;
        silent = true,
        loss_tol = 0.,
        scan_tol = 0.,
        max_iter = 10,
        local_alg = :LD_MMA,
        loss_grad = grad_2p
    )
    @test res2.ret == :MAXEVAL_REACHED

    res3 = get_optimal(
        [4., 1., 1.],
        f_2p;
        silent = true,
        loss_tol = 1e-3,
        scan_tol = 0.,
        local_alg = :LD_MMA,
        scale = [:log, :log, :log],
        loss_grad = grad_2p
    )
    @test isapprox(res3.loss, 5., atol=1e-3)
    @test res3.ret == :FTOL_REACHED
end
