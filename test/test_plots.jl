include("./cases_func.jl")

using ParametersIdentification, Plots

params_plot([3.,4.], 1, f_2p, (1.,6.))

res_f_5p_3im = [params_intervals(
    [3., 0.1, 4, 1.1, 8.],
    i,
    9.,
    f_5p_3im,
    logscale_all = true,
    method = :ONE_PASS
) for i in 1:5]

grid = [params_plot(res_f_5p_3im[i]; max_recursions=5) for i in eachindex(res_f_5p_3im)]

graf = plot(grid[5], label="parameter θ5");
hline!([res_f_5p_3im[5].input.loss_crit], label="confidence level α")
scatter!(grid[5], label="plot grid")
