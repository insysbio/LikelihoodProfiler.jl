using LikelihoodProfiler

include("./cases_func.jl")

#include("../src/params_intervals.jl")
# using ParametersIdentification

res1 = params_intervals(
    [3., 4, 1.1, 10.],
    3,
    9.,
    f_3p_1im,
    logscale_all = true,
    method = :ONE_PASS
)

res2 = params_intervals(
    [3., 4, 1.1, 10.],
    1,
    9.,
    f_3p_1im,
    logscale_all = false,
    method = :D2D_PLE
)
