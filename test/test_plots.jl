using ParametersIdentification

include("./test_func.jl")

params_plot([3.,4.], 1, 8.0, f_2p, (1.,6.))
#params_plot([100.], 1, 8.0, f_1p_ex, (97.,102.))
