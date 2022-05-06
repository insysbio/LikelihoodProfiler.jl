# f_3p_1im_dep(x) = 5.0 + (x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2 # [3., 2., any]

### base example

@testset "default options" begin
    res0 = [get_endpoint(
        [3., 2., 2.1],
        i,
        (x::Vector{Float64}) -> f_3p_1im_dep(x),
        :CICO_ONE_PASS;
        loss_crit = 9.,
        silent = true
    ) for i in 1:3]

    @test isapprox(res0[1].value, 5.0, atol=1e-2)
    @test length(res0[1].profilePoints) > 0
    @test res0[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0[1].direction == :right
    @test isapprox(res0[2].value, 2.0+2.0*sqrt(2.), atol=1e-2)
    @test length(res0[2].profilePoints) > 0
    @test res0[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0[2].direction == :right
    @test length(res0[3].profilePoints) == 0
    @test res0[3].status == :SCAN_BOUND_REACHED
    @test res0[3].direction == :right
end

@testset ":left" begin
    res0 = [get_endpoint(
        [3., 2., 2.1],
        i,
        (x::Vector{Float64}) -> f_3p_1im_dep(x),
        :CICO_ONE_PASS,
        :left;
        loss_crit = 9.,
        silent = true
    ) for i in 1:3]

    @test isapprox(res0[1].value, 1.0, atol=1e-2)
    @test length(res0[1].profilePoints) > 0
    @test res0[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0[1].direction == :left
    @test isapprox(res0[2].value, 2.0-2.0*sqrt(2.), atol=1e-2)
    @test length(res0[2].profilePoints) > 0
    @test res0[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0[2].direction == :left
    @test length(res0[3].profilePoints) == 0
    @test res0[3].status == :SCAN_BOUND_REACHED
    @test res0[3].direction == :left
end

@testset ":log" begin
    res0 = [get_endpoint(
        [3., 2., 2.1],
        i,
        (x::Vector{Float64}) -> f_3p_1im_dep(x),
        :CICO_ONE_PASS,
        :right;
        loss_crit = 9.,
        scale = [:log, :direct, :log],
        silent = true
    ) for i in 1:3]

    @test isapprox(log10(res0[1].value), log10(5.), atol=1e-2)
    @test length(res0[1].profilePoints) > 0
    @test res0[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0[1].direction == :right
    @test isapprox(res0[2].value, 2.0+2.0*sqrt(2.), atol=1e-2)
    @test length(res0[2].profilePoints) > 0
    @test res0[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0[2].direction == :right
    @test length(res0[3].profilePoints) == 0
    @test res0[3].status == :SCAN_BOUND_REACHED
    @test res0[3].direction == :right
end

### derivation free

@testset ":LN_NELDERMEAD, no scale" begin
    res0 = get_endpoint(
        [3., 2., 2.1],
        1,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        local_alg = :LN_NELDERMEAD,
        silent = true
    )
    @test res0.status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0.counter < 1000
    @test isapprox(res0.value, 5.; atol=1e-2)
end

@testset ":LN_NELDERMEAD, no scale, !unpredictable" begin
    res0 = get_endpoint(
        [3., 2., 2.1],
        3,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        local_alg = :LN_NELDERMEAD,
        silent = true
    )
    @test res0.status == :SCAN_BOUND_REACHED
    @test res0.counter < 1000
end

@testset ":LN_NELDERMEAD, log scale" begin
    res0 = get_endpoint(
        [3., 2., 2.1],
        1,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        scale = [:log, :log, :log],
        local_alg = :LN_NELDERMEAD,
        silent = true
    )
    @test res0.status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0.counter < 1000
    @test isapprox(res0.value, 5.; atol=1e-2)
end

@testset ":LN_NELDERMEAD, log scale, !unpredictable" begin
    res0 = get_endpoint(
        [3., 2., 2.1],
        3,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        scale = [:log, :log, :log],
        local_alg = :LN_NELDERMEAD,
        silent = true
    )
    @test res0.status == :SCAN_BOUND_REACHED
    @test res0.counter < 1000
end

### gradient-based

@testset ":LD_SLSQP, no scale, no loss_grad => error" begin
    @test_throws ArgumentError get_endpoint(
        [3., 2., 2.1],
        1,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        local_alg = :LD_SLSQP,
        silent = true
    )
end

@testset ":LD_SLSQP, log scale, no loss_grad => error" begin
    @test_throws ArgumentError get_endpoint(
        [3., 2., 2.1],
        1,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        scale = [:log, :log, :log],
        local_alg = :LD_SLSQP,
        silent = true
    )
end

## autodiff

@testset ":LD_SLSQP, no scale, :AUTODIFF" begin
    res0 = get_endpoint(
        [3., 2., 2.1],
        1,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        local_alg = :LD_SLSQP,
        loss_grad = :AUTODIFF,
        silent = true
    )
    @test res0.status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0.counter < 1000
    @test isapprox(res0.value, 5.; atol=1e-2)
end

@testset ":LD_SLSQP, no scale, :AUTODIFF, !unpredictable" begin
    res0 = get_endpoint(
        [3., 2., 2.1],
        3,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        local_alg = :LD_SLSQP,
        loss_grad = :AUTODIFF,
        silent = true
    )
    @test res0.status == :SCAN_BOUND_REACHED
    @test res0.counter < 1000
end

@testset ":LD_SLSQP, log scale, :AUTODIFF" begin # XXX: scan_tol should be updated
    res0 = get_endpoint(
        [3., 2., 2.1],
        1,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        scale = [:log, :log, :log],
        local_alg = :LD_SLSQP,
        loss_grad = :AUTODIFF,
        scan_tol = 1e-4,
        silent = true
    )
    @test res0.status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0.counter < 1000
    @test isapprox(res0.value, 5.; atol=1e-2)
end

@testset ":LD_SLSQP, log scale, :AUTODIFF, !unpredictable" begin
    res0 = get_endpoint(
        [3., 2., 2.1],
        3,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        scale = [:log, :log, :log],
        local_alg = :LD_SLSQP,
        loss_grad = :AUTODIFF,
        silent = true
    )
    @test res0.status == :SCAN_BOUND_REACHED
    @test res0.counter < 1000
end

## FINITE

@testset ":LD_SLSQP, no scale, :FINITE" begin
    res0 = get_endpoint(
        [3., 2., 2.1],
        1,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        local_alg = :LD_SLSQP,
        loss_grad = :FINITE,
        silent = true
    )
    @test res0.status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0.counter < 1000
    @test isapprox(res0.value, 5.; atol=1e-2)
end

@testset ":LD_SLSQP, no scale, :FINITE, !unpredictable" begin
    res0 = get_endpoint(
        [3., 2., 2.1],
        3,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        local_alg = :LD_SLSQP,
        loss_grad = :FINITE,
        silent = true
    )
    @test res0.status == :SCAN_BOUND_REACHED
    @test res0.counter < 1000
end

@testset ":LD_SLSQP, log scale, :FINITE" begin # XXX: scan_tol should be updated
    res0 = get_endpoint(
        [3., 2., 2.1],
        1,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        scale = [:log, :log, :log],
        local_alg = :LD_SLSQP,
        loss_grad = :FINITE,
        scan_tol = 1e-4,
        silent = true
    )
    @test res0.status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0.counter < 1000
    @test isapprox(log10(res0.value), log10(5.); atol=1e-2)
end

@testset ":LD_SLSQP, log scale, :FINITE, !unpredictable" begin
    res0 = get_endpoint(
        [3., 2., 2.1],
        3,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        scale = [:log, :log, :log],
        local_alg = :LD_SLSQP,
        loss_grad = :FINITE,
        silent = true
    )
    @test res0.status == :SCAN_BOUND_REACHED
    @test res0.counter < 1000
end

## explicit gradient

@testset ":LD_SLSQP, no scale, Function" begin
    res0 = get_endpoint(
        [3., 2., 2.1],
        1,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        local_alg = :LD_SLSQP,
        loss_grad = (x) -> [  # (x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2
            2.0 * (x[1]-3.0) + 2.0 * (x[1]-x[2]-1.0),
            -2.0 * (x[1]-x[2]-1.0),
            0.0
        ], 
        silent = true
    )
    @test res0.status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0.counter < 1000
    @test isapprox(res0.value, 5.; atol=1e-2)
end

@testset ":LD_SLSQP, no scale, Function, !unpredictable" begin
    res0 = get_endpoint(
        [3., 2., 2.1],
        3,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        local_alg = :LD_SLSQP,
        loss_grad = (x) -> [  # (x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2
            2.0 * (x[1]-3.0) + 2.0 * (x[1]-x[2]-1.0),
            -2.0 * (x[1]-x[2]-1.0),
            0.0
        ],
        silent = true
    )
    @test res0.status == :SCAN_BOUND_REACHED
    @test res0.counter < 1000
end

@testset ":LD_SLSQP, log scale, Function" begin
    res0 = get_endpoint(
        [3., 2., 2.1],
        1,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        scale = [:log, :log, :log],
        local_alg = :LD_SLSQP,
        loss_grad = (x) -> [  # (x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2
            2.0 * (x[1]-3.0) + 2.0 * (x[1]-x[2]-1.0),
            -2.0 * (x[1]-x[2]-1.0),
            0.0
        ],
        scan_tol = 1e-4,
        silent = true
    )
    @test res0.status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0.counter < 1000
    @test isapprox(res0.value, 5.; atol=1e-2)
end

@testset ":LD_SLSQP, log scale, Function, !unpredictable" begin
    res0 = get_endpoint(
        [3., 2., 2.1],
        3,            # parameter number
        f_3p_1im_dep, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 9.,
        scale = [:log, :log, :log],
        local_alg = :LD_SLSQP,
        loss_grad = (x) -> [  # (x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2
            2.0 * (x[1]-3.0) + 2.0 * (x[1]-x[2]-1.0),
            -2.0 * (x[1]-x[2]-1.0),
            0.0
        ],
        silent = true
    )
    @test res0.status == :SCAN_BOUND_REACHED
    @test res0.counter < 1000
end