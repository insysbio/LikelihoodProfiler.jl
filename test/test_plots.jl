using Plots
# plotly()

method = :LIN_EXTRAPOL

@testset "plot with no errors" begin
    res0 = [get_interval(
        [3., 2., 2.1],
        i,
        (x::Vector{Float64}) -> f_3p_1im_dep(x),
        method;
        loss_crit = 9.
    ) for i in 1:3]
    p = plot(res0[1])
    @test p isa Plots.Plot
end
