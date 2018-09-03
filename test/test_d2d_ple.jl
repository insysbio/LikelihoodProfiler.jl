using ParametersIdentification

include("./test_func.jl")

params_intervals_d2d([3. ,4.], 1, 8.0, f_2p)
#params_plot([100.], 1, 8.0, f_1p_ex, (97.,102.))
