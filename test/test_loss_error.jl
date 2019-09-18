function err_fun_generate()
    counter::Integer = 0
    fun = (x::Vector{Float64}) -> begin
        counter += 1
        if counter > 5
            throw(ErrorException("Function can be called only 5 times."))
        end
        5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2
    end
end

@testset "loss error in :CICO_ONE_PASS" begin
    err_fun = err_fun_generate()
    res0 = get_endpoint(
        [3., 4.],
        1,
        err_fun,
        :CICO_ONE_PASS;
        loss_crit = 9.
    )
    @test res0.value == nothing
    @test res0.profilePoints == ProfilePoint[]
    @test res0.status == :LOSS_ERROR_STOP
    @test res0.counter < 5
    @test typeof(res0.supreme) == Float64
end

@testset "test profile error" begin
    err_func = err_fun_generate()
    prof = profile(
        [3., 4.],
        1,
        err_func
    )
    res0 = prof(5.)
    @test res0.ret == :FORCED_STOP
    @test res0.counter == 5
end

@testset "loss error in :LIN_EXTRAPOL" begin
    err_fun = err_fun_generate()
    res0 = get_endpoint(
        [3., 4.],
        1,
        err_fun,
        :LIN_EXTRAPOL;
        loss_crit = 9.
    )
    @test res0.value == nothing
    @test length(res0.profilePoints) == 1
    @test res0.status == :LOSS_ERROR_STOP
    @test res0.counter < 5
    @test typeof(res0.supreme) == Float64
    pp = res0.profilePoints[1]
    @test pp.ret == :FORCED_STOP
end

@testset "loss error in :QUADR_EXTRAPOL" begin
    err_fun = err_fun_generate()
    res0 = get_endpoint(
        [3., 4.],
        1,
        err_fun,
        :QUADR_EXTRAPOL;
        loss_crit = 9.
    )
    @test res0.value == nothing
    @test length(res0.profilePoints) == 1
    @test res0.status == :LOSS_ERROR_STOP
    @test res0.counter < 5
    @test typeof(res0.supreme) == Float64
    pp = res0.profilePoints[1]
    @test pp.ret == :FORCED_STOP
end
