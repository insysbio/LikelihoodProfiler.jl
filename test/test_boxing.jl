@testset "boxing" begin
    @testset "1" begin
        box = [
            (-2., 1.1),
            (-2., 1.1),
            (-10., 10.)
        ]

        out = boxing([1.1, 0.0, -1000.], box)

        @test isapprox(out, [1.1, 0.0, -10.])
    end
    @testset "2" begin
        box = [
            (-2., 1.1),
            (-2., 1.1),
            (-Inf, 10.)
        ]

        out = boxing([1.1, 0.0, -1000.], box)

        @test isapprox(out, [1.1, 0.0, -1000.])
    end
end
