using ParametersIdentification, Test

include("./cases_func.jl")

res_f_1p = params_intervals(
    [3.],
    1,
    9.,
    f_1p,
    logscale_all = false,
    method = :ONE_PASS
)

res_f_2p_1im = [params_intervals(
    [3., 1.],
    i,
    9.,
    f_2p_1im,
    logscale_all = false,
    method = :ONE_PASS
) for i in 1:2]

res_f_2p = [params_intervals(
    [3., 4.],
    i,
    9.,
    f_2p,
    logscale_all = false,
    method = :ONE_PASS
) for i in 1:2]

res_f_3p_1im = [params_intervals(
    [3., 4, 1.1],
    i,
    9.,
    f_3p_1im,
    logscale_all = true,
    method = :ONE_PASS
) for i in 1:3]

res_f_4p_2im = [params_intervals(
    [3., 4, 1.1, 10.],
    i,
    9.,
    f_4p_2im,
    logscale_all = true,
    method = :ONE_PASS
) for i in 1:4]

res_f_4p_3im = [params_intervals(
    [3., 4, 1.1, 10.],
    i,
    9.,
    f_4p_3im,
    logscale_all = true,
    method = :ONE_PASS
) for i in 1:4]

# tests
@testset "f_1p" begin
    @test all(@. isapprox(res_f_1p.intervals, [1.0, 5.0], atol=0.001))
    @test res_f_1p.ret_codes == [:FTOL_REACHED, :FTOL_REACHED]
end

@testset "f_2p_1im" begin
    @test all(@. isapprox(res_f_2p_1im[1].intervals, [1.0, 5.0], atol=0.001))
    @test res_f_2p_1im[1].ret_codes == [:FTOL_REACHED, :FTOL_REACHED]
    @test res_f_2p_1im[2].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
end

@testset "f_2p" begin
    @test all(@. isapprox(res_f_2p[1].intervals, [1.0, 5.0], atol=0.001))
    @test res_f_2p[1].ret_codes == [:FTOL_REACHED, :FTOL_REACHED]
    @test all(@. isapprox(res_f_2p[2].intervals, [2.0, 6.0], atol=0.001))
    @test res_f_2p[2].ret_codes == [:FTOL_REACHED, :FTOL_REACHED]
end

@testset "f_3p_1im" begin
    @test all(@. isapprox(res_f_3p_1im[1].intervals, [1.0, 5.0], atol=0.001))
    @test res_f_3p_1im[1].ret_codes == [:FTOL_REACHED, :FTOL_REACHED]
    @test res_f_3p_1im[2].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
    @test res_f_3p_1im[3].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
end

@testset "f_4p_2im" begin
    @test all(@. isapprox(res_f_4p_2im[1].intervals, [1.0, 5.0], atol=0.001))
    @test res_f_4p_2im[1].ret_codes == [:FTOL_REACHED, :FTOL_REACHED]
    @test all(@. isapprox(res_f_4p_2im[2].intervals, [2.0, 6.0], atol=0.001))
    @test res_f_4p_2im[2].ret_codes == [:FTOL_REACHED, :FTOL_REACHED]
    @test res_f_4p_2im[3].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
    @test res_f_4p_2im[4].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
end

@testset "f_4p_3im" begin
    @test all(@. isapprox(res_f_4p_3im[1].intervals, [1.0, 5.0], atol=0.001))
    @test res_f_4p_3im[1].ret_codes == [:FTOL_REACHED, :FTOL_REACHED]
    @test res_f_4p_3im[2].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
    @test res_f_4p_3im[3].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
    @test res_f_4p_3im[4].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
end
