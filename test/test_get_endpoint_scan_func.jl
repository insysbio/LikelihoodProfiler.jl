### derivation free

@testset ":LN_NELDERMEAD, no scale" begin
    res0 = get_endpoint(
        [2.,2.,2.],
        (x) -> x[1]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        local_alg = :LN_NELDERMEAD
    )
    @test res0.status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0.counter < 1000
    @test isapprox(res0.value, 169.; atol=1e-2)
end

@testset ":LN_NELDERMEAD, no scale, !unpredictable" begin
    res0 = get_endpoint(
        [2.,2.,2.],
        (x) -> x[3]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        scan_bound = 200.,
        local_alg = :LN_NELDERMEAD
    )
    @test res0.status == :SCAN_BOUND_REACHED
    @test res0.counter < 1000
end

@testset ":LN_NELDERMEAD, log scale" begin
    res0 = get_endpoint(
        [2.,2.,2.],
        (x) -> x[1]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        scale = [:log, :log, :log],
        local_alg = :LN_NELDERMEAD
    )
    @test res0.status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0.counter < 1000
    @test isapprox(res0.value, 169.; atol=1e-2)
end

@testset ":LN_NELDERMEAD, log scale, !unpredictable" begin
    res0 = get_endpoint(
        [2.,2.,2.],
        (x) -> x[3]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        scale = [:log, :log, :log],
        scan_bound = 200.,
        local_alg = :LN_NELDERMEAD
    )
    @test res0.status == :SCAN_BOUND_REACHED
    @test res0.counter < 1000
end

### gradient-based

@testset ":LD_SLSQP, no scale, no scan_grad => error" begin
    @test_throws ArgumentError get_endpoint(
        [2.,2.,2.],
        (x) -> x[1]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        local_alg = :LD_SLSQP
    )
end

@testset ":LD_SLSQP, log scale, no scan_grad => error" begin
    @test_throws ArgumentError get_endpoint(
        [2.,2.,2.],
        (x) -> x[1]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        scale = [:log, :log, :log],
        local_alg = :LD_SLSQP
    )
end

@testset ":LD_SLSQP, no scale, :AUTODIFF, :AUTODIFF" begin
    res0 = get_endpoint(
        [2.,2.,2.],
        (x) -> x[1]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        local_alg = :LD_SLSQP,
        scan_grad = :AUTODIFF,
        loss_grad = :AUTODIFF
    )
    @test res0.status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0.counter < 1000
    @test isapprox(res0.value, 169.; atol=1e-2)
end

@testset ":LD_SLSQP, no scale, :AUTODIFF, :AUTODIFF, !unpredictable" begin
    res0 = get_endpoint(
        [2.,2.,2.],
        (x) -> x[3]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        scan_bound = 200.,
        local_alg = :LD_SLSQP,
        scan_grad = :AUTODIFF,
        loss_grad = :AUTODIFF
    )
    @test res0.status == :SCAN_BOUND_REACHED
    @test res0.counter < 1000
end

@testset ":LD_SLSQP, log scale, :AUTODIFF, :AUTODIFF" begin
    res0 = get_endpoint(
        [2.,2.,2.],
        (x) -> x[1]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        scale = [:log, :log, :log],
        local_alg = :LD_SLSQP,
        scan_grad = :AUTODIFF,
        loss_grad = :AUTODIFF
    )
    @test res0.status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0.counter < 1000
    @test isapprox(res0.value, 169.; atol=1e-2)
end

@testset ":LD_SLSQP, log scale, :AUTODIFF, :AUTODIFF, !unpredictable" begin
    res0 = get_endpoint(
        [2.,2.,2.],
        (x) -> x[3]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        scan_bound = 200.,
        scale = [:log, :log, :log],
        local_alg = :LD_SLSQP,
        scan_grad = :AUTODIFF,
        loss_grad = :AUTODIFF
    )
    @test res0.status == :SCAN_BOUND_REACHED
    @test res0.counter < 1000
end

@testset ":LD_SLSQP, no scale, :FINITE, :FINITE" begin
    res0 = get_endpoint(
        [2.,2.,2.],
        (x) -> x[1]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        local_alg = :LD_SLSQP,
        scan_grad = :FINITE,
        loss_grad = :FINITE
    )
    @test res0.status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0.counter < 1000
    @test isapprox(res0.value, 169.; atol=1e-2)
end

@testset ":LD_SLSQP, no scale, :FINITE, :FINITE, !unpredictable" begin
    res0 = get_endpoint(
        [2.,2.,2.],
        (x) -> x[3]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        scan_bound = 200.,
        local_alg = :LD_SLSQP,
        scan_grad = :FINITE,
        loss_grad = :FINITE
    )
    @test res0.status == :SCAN_BOUND_REACHED
    @test res0.counter < 1000
end

@testset ":LD_SLSQP, log scale, :FINITE, :FINITE" begin
    res0 = get_endpoint(
        [2.,2.,2.],
        (x) -> x[1]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        scale = [:log, :log, :log],
        local_alg = :LD_SLSQP,
        scan_grad = :FINITE,
        loss_grad = :FINITE
    )
    @test res0.status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0.counter < 1000
    @test isapprox(res0.value, 169.; atol=1e-2)
end

@testset ":LD_SLSQP, log scale, :FINITE, :FINITE, !unpredictable" begin
    res0 = get_endpoint(
        [2.,2.,2.],
        (x) -> x[3]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        scan_bound = 200.,
        scale = [:log, :log, :log],
        local_alg = :LD_SLSQP,
        scan_grad = :FINITE,
        loss_grad = :FINITE
    )
    @test res0.status == :SCAN_BOUND_REACHED
    @test res0.counter < 1000
end

@testset ":LD_SLSQP, no scale, Function, Function" begin
    res0 = get_endpoint(
        [2.,2.,2.],
        (x) -> x[1]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        local_alg = :LD_SLSQP,
        scan_grad = (x) -> [2.0 * x[1], 0., 0.],
        loss_grad = (x) -> [2.0 * (x[1] - 5.), 2.0 * (x[2] - 6.), 0.]
    )
    @test res0.status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0.counter < 1000
    @test isapprox(res0.value, 169.; atol=1e-2)
end

@testset ":LD_SLSQP, no scale, Function, Function, !unpredictable" begin
    res0 = get_endpoint(
        [2.,2.,2.],
        (x) -> x[3]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        scan_bound = 200.,
        local_alg = :LD_SLSQP,
        scan_grad = (x) -> [0., 0., 2.0 * x[3]],
        loss_grad = (x) -> [2.0 * (x[1] - 5.), 2.0 * (x[2] - 6.), 0.]
    )
    @test res0.status == :SCAN_BOUND_REACHED
    @test res0.counter < 1000
end

@testset ":LD_SLSQP, log scale, Function, Function" begin
    res0 = get_endpoint(
        [2.,2.,2.],
        (x) -> x[1]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        scale = [:log, :log, :log],
        local_alg = :LD_SLSQP,
        scan_grad = (x) -> [2.0 * x[1], 0., 0.],
        loss_grad = (x) -> [2.0 * (x[1] - 5.), 2.0 * (x[2] - 6.), 0.]
    )
    @test res0.status == :BORDER_FOUND_BY_SCAN_TOL
    @test res0.counter < 1000
    @test isapprox(res0.value, 169.; atol=1e-2)
end

@testset ":LD_SLSQP, log scale, Function, Function, !unpredictable" begin
    res0 = get_endpoint(
        [2.,2.,2.],
        (x) -> x[3]^2, # scan func
        (x) -> (x[1]-5.)^2 + (x[2] - 6.)^2, # loss_func
        :CICO_ONE_PASS,
        loss_crit = 64.,
        scan_bound = 200.,
        scale = [:log, :log, :log],
        local_alg = :LD_SLSQP,
        scan_grad = (x) -> [0., 0., 2.0 * x[3]],
        loss_grad = (x) -> [2.0 * (x[1] - 5.), 2.0 * (x[2] - 6.), 0.]
    )
    @test res0.status == :SCAN_BOUND_REACHED
    @test res0.counter < 1000
end
