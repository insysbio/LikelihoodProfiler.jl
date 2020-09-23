#=

ADNLPModel(f, x0)
ADNLPModel(f, x0, lvar, uvar)
ADNLPModel(f, x0, c, lcon, ucon)
ADNLPModel(f, x0, lvar, uvar, c, lcon, ucon)

ADNLPModel is an AbstractNLPModel using ForwardDiff to compute the derivatives. The problem is defined as

 min  f(x)
s.to  lcon ≤ c(x) ≤ ucon
      lvar ≤   x  ≤ uvar.

The following keyword arguments are available to all constructors:
  •    name: The name of the model (default: "Generic")
The following keyword arguments are available to the constructors for constrained problems:
  •    lin: An array of indexes of the linear constraints (default: Int[])
  •    y0: An inital estimate to the Lagrangian multipliers (default: zeros)
 =#

using NLPModels, Percival

f0(x) = x[1]
f1(x) = -x[1]^2
f2(x) = log10(x[1]^2)
x01 = [3.0]
x02 = [3.0, 2.0, 1.0]
lvar1 = [1e-9]
uvar1 = [1e9]
lvar2 = [1e-9,1e-9,1e-9]
uvar2 = [1e9,1e9,1e9]
lcon = [-Inf]
ucon = [0.]
c1(x) = [(x[1]-3.0)^2-4.0]
c2(x) =[(x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2 - 4.0]

m0 = ADNLPModel(f0, x01, lvar1, uvar1, c1, lcon, ucon)
m1 = ADNLPModel(f1, x0, lvar, uvar, c1, lcon, ucon)
m2 = ADNLPModel(f2, x0, lvar, uvar, c1, lcon, ucon)

output1 = percival(m0,subsolver=:LN_NELDERMEAD)
output2 = percival(m2, subsolver=:nelder_mead)

println(output1)
println(output2)

using Test

@testset "f_2p_1im" begin
    ret = [get_interval_percival(
       [3.0; 1.0],
       i,
       f_2p_1im,
       9.0,
    ) for i in 1:2]
    @test ret[1][1].solution[1] ≈ 1.0 rtol=1e-6
    @test ret[1][2].solution[1] ≈ 5.0 rtol=1e-6
    @test ret[2][1].solution[2] ≈ 1e-9 rtol=1e-6
    @test ret[2][2].solution[2] ≈ 1e9 rtol=1e-6
end

@testset "f_2p" begin
    ret = [get_interval_percival(
       [3.0; 4.0],
       i,
       f_2p,
       9.0,
    ) for i in 1:2]
    @test ret[1][1].solution[1] ≈ 1.0 rtol=1e-6
    @test ret[1][2].solution[1] ≈ 5.0 rtol=1e-6
    @test ret[2][1].solution[2] ≈ 2.0 rtol=1e-6
    @test ret[2][2].solution[2] ≈ 6.0 rtol=1e-6
end

@testset "f_3p_1im" begin
    ret = [get_interval_percival(
       [3.0; 1.0; 1.0],
       i,
       f_3p_1im,
       9.0,
    ) for i in 1:3]
    @test ret[1][1].solution[1] ≈ 1.0 rtol=1e-6
    @test ret[1][2].solution[1] ≈ 5.0 rtol=1e-6
    @test ret[2][1].solution[2] ≈ 1e-9 rtol=1e-6
    @test ret[2][2].solution[2] ≈ 1e9 rtol=1e-6
    @test ret[3][1].solution[3] ≈ 1e-9 rtol=1e-6
    @test ret[3][2].solution[3] ≈ 1e9 rtol=1e-6
end

@testset "f_3p_1im_dep" begin
    ret = jump_identify(
        f_3p_1im_dep,
        [3.0, 1.0, 1.0],
        scan_bounds=(lb,ub),
        solver=solver
    )
    @test ret[1].ci[1] ≈ 1.0 atol=tol
    @test ret[1].ci[2] ≈ 5.0 atol=tol
    @test ret[2].ci[1] ≈ lb atol=tol
    @test ret[2].ci[2] ≈ 2.0+2.0*sqrt(2.) atol=tol
    @test ret[3].ci[1] ≈ lb atol=tol
    @test ret[3].ci[2] ≈ ub rtol=tol
end

@testset "f_4p_2im" begin
    ret = jump_identify(
        f_4p_2im,
        [3.0, 4.0, 1.0, 1.0],
        scan_bounds=(lb,ub),
        solver=solver
    )
    @test ret[1].ci[1] ≈ 1.0 atol=tol
    @test ret[1].ci[2] ≈ 5.0 atol=tol
    @test ret[2].ci[1] ≈ 2.0 atol=tol
    @test ret[2].ci[2] ≈ 6.0 atol=tol
    @test ret[3].ci[1] ≈ lb atol=tol
    @test ret[3].ci[2] ≈ ub rtol=tol
    @test ret[4].ci[1] ≈ lb atol=tol
    @test ret[4].ci[2] ≈ ub rtol=tol
end

@testset "f_4p_3im" begin
    ret = jump_identify(
        f_4p_3im,
        [3.0, 4.0, 1.0, 1.0],
        scan_bounds=(lb,ub),
        solver=solver
    )
    @test ret[1].ci[1] ≈ 1.0 atol=tol
    @test ret[1].ci[2] ≈ 5.0 atol=tol
    @test ret[2].ci[1] ≈ lb atol=tol
    @test ret[2].ci[2] ≈ ub rtol=tol
    @test ret[3].ci[1] ≈ lb atol=tol
    @test ret[3].ci[2] ≈ ub rtol=tol #failed NUMERICAL_ERROR
    @test ret[4].ci[1] ≈ lb atol=tol
    @test ret[4].ci[2] ≈ ub rtol=tol
end

#=
@testset "f_1p_ex" begin
    ret = jump_identify(
        f_1p_ex,
        [1e-8]
    )
    @test ret[1][1] ≈ 1e-8-2.0 atol=tol
    @test ret[1][2] ≈ 1e-8+2.0 atol=tol
end
=#

@testset "f_5p_3im" begin
    ret = jump_identify(
        f_5p_3im,
        [3.0, 0.0, 1.0, 1.0, 1.0],
        scan_bounds=(lb,ub),
        solver=solver
    )
    @test ret[1].ci[1] ≈ 1.0 atol=tol
    @test ret[1].ci[2] ≈ 5.0 atol=tol
    @test ret[2].ci[1] ≈ lb atol=tol
    @test ret[2].ci[2] ≈ log(3.) atol=tol
    @test ret[3].ci[1] ≈ lb atol=tol #failed NUMERICAL_ERROR
    @test ret[3].ci[2] ≈ ub rtol=tol
    @test ret[4].ci[1] ≈ lb atol=tol
    @test ret[4].ci[2] ≈ ub rtol=tol #LOCALLY_INFEASIBLE
    @test ret[5].ci[1] ≈ lb atol=tol
    @test ret[5].ci[2] ≈ ub rtol=tol #failed NUMERICAL_ERROR
end

@testset "f_3p_im" begin
    ret = jump_identify(
        f_3p_im,
        [3.0, 1.0, 1.0],
        scan_bounds=(lb,ub),
        solver=solver
    )
    @test ret[1].ci[1] ≈ 1.0 atol=tol
    @test ret[1].ci[2] ≈ 5.0 atol=tol
    @test ret[2].ci[1] ≈ lb atol=tol
    @test ret[2].ci[2] ≈ log(3.) atol=tol
    @test ret[3].ci[1] ≈ lb atol=tol
    @test ret[3].ci[2] ≈ ub rtol=tol
end
