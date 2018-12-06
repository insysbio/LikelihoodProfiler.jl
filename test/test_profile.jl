@testset "profile" begin
    prof = profile(
        [3., 2., 2.1],
        1,
        f_3p_1im_dep;
        skip_optim = false
    )

    x =  0.:0.2:10
    y = [prof(x[i]) for i in 1:length(x)]
    
    @test isapprox(y[5].value, 0.8)
    @test isapprox(y[5].loss, 9.8402, atol=1e-3)
    @test y[5].ret == :FTOL_REACHED
    @test y[5].counter > 1
    @test length(y[5].params) == 3
end
