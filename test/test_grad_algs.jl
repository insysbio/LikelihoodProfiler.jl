### Gradient-based optimization tests
# The following NLopt gradient -based algorithms are compared
grad_algorithms = [
    :LD_MMA, # Method of Moving Asymptotes
    :LD_SLSQP, # Sequential Least-Squares Quadratic Programming
    :LD_LBFGS, # Low-storage BFGS
    :LD_TNEWTON_PRECOND_RESTART, # Preconditioned truncated Newton
    :LD_TNEWTON_PRECOND, # Same without restarting
    :LD_TNEWTON_RESTART, # Same without preconditioning
    :LD_TNEWTON, # Same without restarting or preconditioning
    :LD_VAR2, # Shifted limited-memory variable-metric (rank 2)
    :LD_VAR1  # Shifted limited-memory variable-metric (rank 1)
]

grad_res = Dict{Symbol, Vector{ParamInterval}}()

for alg in grad_algorithms
    grad_res[alg] = [get_interval(
        [3., 2., 2.1],
        i,
        f_3p_1im_dep,
        :CICO_ONE_PASS,
        local_alg = alg,
        loss_crit = 9.
    ) for i in 1:3]
end

@testset "LD_MMA" begin
    alg = :LD_MMA
    @test isapprox(grad_res[alg][1].result[1].value, 1.0, atol=1e-2)
    @test isapprox(grad_res[alg][1].result[2].value, 5.0, atol=1e-2)
    @test length(grad_res[alg][1].result[1].profilePoints) > 0
    @test length(grad_res[alg][1].result[2].profilePoints) > 0
    @test grad_res[alg][1].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[1].direction == :left
    @test grad_res[alg][1].result[2].direction == :right
    @test isapprox(grad_res[alg][2].result[1].value, 2.0-2.0*sqrt(2.), atol=1e-2)
    @test isapprox(grad_res[alg][2].result[2].value, 2.0+2.0*sqrt(2.), atol=1e-2)
    @test grad_res[alg][2].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][2].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][3].result[1].status == :SCAN_BOUND_REACHED
    @test grad_res[alg][3].result[2].status == :SCAN_BOUND_REACHED
end

@testset "LD_SLSQP" begin
    alg = :LD_SLSQP
    @test isapprox(grad_res[alg][1].result[1].value, 1.0, atol=1e-2)
    @test isapprox(grad_res[alg][1].result[2].value, 5.0, atol=1e-2)
    @test length(grad_res[alg][1].result[1].profilePoints) > 0
    @test length(grad_res[alg][1].result[2].profilePoints) > 0
    @test grad_res[alg][1].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[1].direction == :left
    @test grad_res[alg][1].result[2].direction == :right
    @test isapprox(grad_res[alg][2].result[1].value, 2.0-2.0*sqrt(2.), atol=1e-2)
    @test isapprox(grad_res[alg][2].result[2].value, 2.0+2.0*sqrt(2.), atol=1e-2)
    @test grad_res[alg][2].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][2].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][3].result[1].status == :SCAN_BOUND_REACHED
    @test grad_res[alg][3].result[2].status == :SCAN_BOUND_REACHED
end

@testset "LD_LBFGS" begin
    alg = :LD_LBFGS
    @test isapprox(grad_res[alg][1].result[1].value, 1.0, atol=1e-2)
    @test isapprox(grad_res[alg][1].result[2].value, 5.0, atol=1e-2)
    @test length(grad_res[alg][1].result[1].profilePoints) > 0
    @test length(grad_res[alg][1].result[2].profilePoints) > 0
    @test grad_res[alg][1].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[1].direction == :left
    @test grad_res[alg][1].result[2].direction == :right
    @test isapprox(grad_res[alg][2].result[1].value, 2.0-2.0*sqrt(2.), atol=1e-2)
    @test isapprox(grad_res[alg][2].result[2].value, 2.0+2.0*sqrt(2.), atol=1e-2)
    @test grad_res[alg][2].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][2].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][3].result[1].status == :SCAN_BOUND_REACHED
    @test grad_res[alg][3].result[2].status == :SCAN_BOUND_REACHED
end

@testset "LD_TNEWTON_PRECOND_RESTART" begin
    alg = :LD_TNEWTON_PRECOND_RESTART
    @test isapprox(grad_res[alg][1].result[1].value, 1.0, atol=1e-2)
    @test isapprox(grad_res[alg][1].result[2].value, 5.0, atol=1e-2)
    @test length(grad_res[alg][1].result[1].profilePoints) > 0
    @test length(grad_res[alg][1].result[2].profilePoints) > 0
    @test grad_res[alg][1].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[1].direction == :left
    @test grad_res[alg][1].result[2].direction == :right
    @test isapprox(grad_res[alg][2].result[1].value, 2.0-2.0*sqrt(2.), atol=1e-2)
    @test isapprox(grad_res[alg][2].result[2].value, 2.0+2.0*sqrt(2.), atol=1e-2)
    @test grad_res[alg][2].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][2].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][3].result[1].status == :SCAN_BOUND_REACHED
    @test grad_res[alg][3].result[2].status == :SCAN_BOUND_REACHED
end

@testset "LD_TNEWTON_PRECOND" begin
    alg = :LD_TNEWTON_PRECOND
    @test isapprox(grad_res[alg][1].result[1].value, 1.0, atol=1e-2)
    @test isapprox(grad_res[alg][1].result[2].value, 5.0, atol=1e-2)
    @test length(grad_res[alg][1].result[1].profilePoints) > 0
    @test length(grad_res[alg][1].result[2].profilePoints) > 0
    @test grad_res[alg][1].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[1].direction == :left
    @test grad_res[alg][1].result[2].direction == :right
    @test isapprox(grad_res[alg][2].result[1].value, 2.0-2.0*sqrt(2.), atol=1e-2)
    @test isapprox(grad_res[alg][2].result[2].value, 2.0+2.0*sqrt(2.), atol=1e-2)
    @test grad_res[alg][2].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][2].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][3].result[1].status == :SCAN_BOUND_REACHED
    @test grad_res[alg][3].result[2].status == :SCAN_BOUND_REACHED
end

@testset "LD_TNEWTON_RESTART" begin
    alg = :LD_TNEWTON_RESTART
    @test isapprox(grad_res[alg][1].result[1].value, 1.0, atol=1e-2)
    @test isapprox(grad_res[alg][1].result[2].value, 5.0, atol=1e-2)
    @test length(grad_res[alg][1].result[1].profilePoints) > 0
    @test length(grad_res[alg][1].result[2].profilePoints) > 0
    @test grad_res[alg][1].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[1].direction == :left
    @test grad_res[alg][1].result[2].direction == :right
    @test isapprox(grad_res[alg][2].result[1].value, 2.0-2.0*sqrt(2.), atol=1e-2)
    @test isapprox(grad_res[alg][2].result[2].value, 2.0+2.0*sqrt(2.), atol=1e-2)
    @test grad_res[alg][2].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][2].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test_broken grad_res[alg][3].result[1].status == :SCAN_BOUND_REACHED
    @test_broken grad_res[alg][3].result[2].status == :SCAN_BOUND_REACHED
end

@testset "LD_TNEWTON" begin
    alg = :LD_TNEWTON
    @test isapprox(grad_res[alg][1].result[1].value, 1.0, atol=1e-2)
    @test isapprox(grad_res[alg][1].result[2].value, 5.0, atol=1e-2)
    @test length(grad_res[alg][1].result[1].profilePoints) > 0
    @test length(grad_res[alg][1].result[2].profilePoints) > 0
    @test grad_res[alg][1].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[1].direction == :left
    @test grad_res[alg][1].result[2].direction == :right
    @test isapprox(grad_res[alg][2].result[1].value, 2.0-2.0*sqrt(2.), atol=1e-2)
    @test isapprox(grad_res[alg][2].result[2].value, 2.0+2.0*sqrt(2.), atol=1e-2)
    @test grad_res[alg][2].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][2].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][3].result[1].status == :SCAN_BOUND_REACHED
    @test grad_res[alg][3].result[2].status == :SCAN_BOUND_REACHED
end

@testset "LD_VAR2" begin
    alg = :LD_VAR2
    @test isapprox(grad_res[alg][1].result[1].value, 1.0, atol=1e-2)
    @test isapprox(grad_res[alg][1].result[2].value, 5.0, atol=1e-2)
    @test length(grad_res[alg][1].result[1].profilePoints) > 0
    @test length(grad_res[alg][1].result[2].profilePoints) > 0
    @test grad_res[alg][1].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[1].direction == :left
    @test grad_res[alg][1].result[2].direction == :right
    @test isapprox(grad_res[alg][2].result[1].value, 2.0-2.0*sqrt(2.), atol=1e-2)
    @test isapprox(grad_res[alg][2].result[2].value, 2.0+2.0*sqrt(2.), atol=1e-2)
    @test grad_res[alg][2].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][2].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test_broken grad_res[alg][3].result[1].status == :SCAN_BOUND_REACHED
    @test_broken grad_res[alg][3].result[2].status == :SCAN_BOUND_REACHED
end

@testset "LD_VAR1" begin
    alg = :LD_VAR1
    @test isapprox(grad_res[alg][1].result[1].value, 1.0, atol=1e-2)
    @test isapprox(grad_res[alg][1].result[2].value, 5.0, atol=1e-2)
    @test length(grad_res[alg][1].result[1].profilePoints) > 0
    @test length(grad_res[alg][1].result[2].profilePoints) > 0
    @test grad_res[alg][1].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][1].result[1].direction == :left
    @test grad_res[alg][1].result[2].direction == :right
    @test isapprox(grad_res[alg][2].result[1].value, 2.0-2.0*sqrt(2.), atol=1e-2)
    @test isapprox(grad_res[alg][2].result[2].value, 2.0+2.0*sqrt(2.), atol=1e-2)
    @test grad_res[alg][2].result[1].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][2].result[2].status == :BORDER_FOUND_BY_SCAN_TOL
    @test grad_res[alg][3].result[1].status == :SCAN_BOUND_REACHED
    @test grad_res[alg][3].result[2].status == :SCAN_BOUND_REACHED
end
