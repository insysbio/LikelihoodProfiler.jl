include("./cases_func.jl")

include("../src/params_intervals.jl")
# using ParametersIdentification

res1 = params_intervals(
    [3., 4, 1.1, 10.],
    1,
    9.,
    f_3p_1im,
    logscale_all = true
)
