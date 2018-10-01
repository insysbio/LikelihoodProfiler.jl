
include("./cases_func.jl")

res_f_1p = params_intervals(
    [3.],
    1,
    9.,
    f_1p,
    logscale_all = false,
    method = :D2D_PLE
)

res_f_2p_1im = [params_intervals(
    [3., 1.],
    i,
    9.,
    f_2p_1im,
    logscale_all = false,
    method = :D2D_PLE
) for i in 1:2]

res_f_2p = [params_intervals(
    [3., 4.],
    i,
    9.,
    f_2p,
    logscale_all = false,
    method = :D2D_PLE
) for i in 1:2]

res_f_3p_1im = [params_intervals(
    [3., 4., 1.1],
    i,
    9.,
    f_3p_1im,
    #max_iter = 1000000,
    logscale_all = false,
    method = :D2D_PLE
) for i in 1:3]

res_f_4p_2im = [params_intervals(
    [3., 4, 1.1, 7.],
    i,
    9.,
    f_4p_2im,
    logscale_all = false,
    method = :D2D_PLE
) for i in 1:4]

res_f_4p_3im = [params_intervals(
    [3., 4, 1.1, 7.],
    i,
    9.,
    f_4p_3im,
    logscale_all = false,
    method = :D2D_PLE
) for i in 1:4]

res_f_5p_3im = [params_intervals(
    [3., 0.1, 4, 1.1, 8.],
    i,
    9.,
    f_5p_3im,
    logscale_all = false,
    method = :D2D_PLE
) for i in 1:5]

res_f_3p_im = [params_intervals(
    [3., 0.1, 8.],
    i,
    9.,
    f_3p_im,
    scan_bound = [0.,10^3],
    logscale_all = false,
    method = :D2D_PLE
) for i in 1:3]
# tests
@testset "f_1p" begin
    @test all(@. isapprox(res_f_1p.interval, [1.0, 5.0], atol=res_f_1p.input.ptol))
    @test all(@. isapprox(res_f_1p.loss_final, [res_f_1p.input.loss_crit, res_f_1p.input.loss_crit], atol=res_f_1p.input.losstol))
    @test res_f_1p.ret_codes == [:FTOL_REACHED, :FTOL_REACHED]
end

@testset "f_2p_1im" begin
    @test all(@. isapprox(res_f_2p_1im[1].interval, [1.0, 5.0], atol=res_f_2p_1im[1].input.ptol))
    @test all(@. isapprox(res_f_2p_1im[1].loss_final, [res_f_2p_1im[1].input.loss_crit, res_f_2p_1im[1].input.loss_crit],
    atol=res_f_2p_1im[1].input.losstol))
    @test res_f_2p_1im[1].ret_codes == [:FTOL_REACHED, :FTOL_REACHED]
    @test res_f_2p_1im[2].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
end

@testset "f_2p" begin
    @test all(@. isapprox(res_f_2p[1].interval, [1.0, 5.0], atol=res_f_2p[1].input.ptol))
    @test all(@. isapprox(res_f_2p[1].loss_final, [res_f_2p[1].input.loss_crit, res_f_2p[1].input.loss_crit],
    atol=res_f_2p[1].input.losstol))
    @test res_f_2p[1].ret_codes == [:FTOL_REACHED, :FTOL_REACHED]
    @test all(@. isapprox(res_f_2p[2].interval, [2.0, 6.0], atol=res_f_2p[2].input.ptol))
    @test all(@. isapprox(res_f_2p[2].loss_final, [res_f_2p[2].input.loss_crit, res_f_2p[2].input.loss_crit],
    atol=res_f_2p[2].input.losstol))
    @test res_f_2p[2].ret_codes == [:FTOL_REACHED, :FTOL_REACHED]
end

@testset "f_3p_1im" begin
    @test all(@. isapprox(res_f_3p_1im[1].interval, [1.0, 5.0], atol=res_f_3p_1im[1].input.ptol))
    @test all(@. isapprox(res_f_3p_1im[1].loss_final, [res_f_3p_1im[1].input.loss_crit, res_f_3p_1im[1].input.loss_crit],
    atol=res_f_3p_1im[1].input.losstol))
    @test res_f_3p_1im[1].ret_codes == [:FTOL_REACHED, :FTOL_REACHED]
    @test res_f_3p_1im[2].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
    @test res_f_3p_1im[3].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
end

@testset "f_4p_2im" begin
    @test all(@. isapprox(res_f_4p_2im[1].interval, [1.0, 5.0], atol=res_f_4p_2im[1].input.ptol))
    @test all(@. isapprox(res_f_4p_2im[1].loss_final, [res_f_4p_2im[1].input.loss_crit, res_f_4p_2im[1].input.loss_crit],
    atol=res_f_4p_2im[1].input.losstol))
    @test res_f_4p_2im[1].ret_codes == [:FTOL_REACHED, :FTOL_REACHED]
    @test all(@. isapprox(res_f_4p_2im[2].interval, [2.0, 6.0], atol=res_f_4p_2im[2].input.ptol))
    @test all(@. isapprox(res_f_4p_2im[2].loss_final, [res_f_4p_2im[2].input.loss_crit, res_f_4p_2im[2].input.loss_crit],
    atol=res_f_4p_2im[2].input.losstol))
    @test res_f_4p_2im[2].ret_codes == [:FTOL_REACHED, :FTOL_REACHED]
    @test res_f_4p_2im[3].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
    @test res_f_4p_2im[4].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
end

@testset "f_4p_3im" begin
    @test all(@. isapprox(res_f_4p_3im[1].interval, [1.0, 5.0], atol=res_f_4p_3im[1].input.ptol))
    @test all(@. isapprox(res_f_4p_3im[1].loss_final, [res_f_4p_3im[1].input.loss_crit, res_f_4p_3im[1].input.loss_crit],
    atol=res_f_4p_3im[1].input.losstol))
    @test res_f_4p_3im[1].ret_codes == [:FTOL_REACHED, :FTOL_REACHED]
    @test res_f_4p_3im[2].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
    @test res_f_4p_3im[3].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
    @test res_f_4p_3im[4].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
end

@testset "f_5p_3im" begin
    @test all(@. isapprox(log10.(res_f_5p_3im[1].interval), log10.([1.0, 5.0]), atol=res_f_5p_3im[1].input.ptol))
    @test all(@. isapprox(res_f_5p_3im[1].loss_final, [res_f_5p_3im[1].input.loss_crit, res_f_5p_3im[1].input.loss_crit],
    atol=res_f_5p_3im[1].input.losstol))
    @test res_f_5p_3im[1].ret_codes == [:FTOL_REACHED, :FTOL_REACHED]
    @test isapprox(log10(res_f_5p_3im[2].interval[2]), log10(1.1), atol=res_f_5p_3im[2].input.ptol)
    @test res_f_5p_3im[2].ret_codes == [:BOUNDS_REACHED, :FTOL_REACHED]
    @test res_f_5p_3im[3].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
    @test res_f_5p_3im[4].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
    @test res_f_5p_3im[5].ret_codes == [:BOUNDS_REACHED, :BOUNDS_REACHED]
end
