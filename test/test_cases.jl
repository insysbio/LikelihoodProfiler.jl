include("./cases_func.jl")

using ParametersIdentification

res1 = params_intervals(
    [3., 4, 1.1, 10.],
    3,
    9.,
    f_3p_1im,
    #tol_glob = 1e-3,
    tol_loc = 1e-3,
    max_iter = 100000,
    bounds_params = fill([0., Inf], 4),
    logscale = fill(true, 4),
    bounds_id = [1e-9, 1e9]
)
