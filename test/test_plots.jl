using LikelihoodProfiler, Plots
plotly()

include("./cases_func.jl")

method = :LIN_EXTRAPOL#:CICO_ONE_PASS#

res0 = [get_interval(
    [3., 2., 2.1],
    i,
    (x::Vector{Float64}) -> f_3p_1im_dep(x),
    method;
    loss_crit = 9.
) for i in 1:3]
