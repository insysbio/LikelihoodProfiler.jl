include("./cases_func.jl")

using LikelihoodProfiler, Plots

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

graf = plot(grid[5], label="parameter θ₅");
hline!([res_f_5p_3im[5].input.loss_crit], label="confidence level α")
scatter!(grid[5], label="plot grid")


grid_dream = [params_plot(r; max_recursions=2) for r in result]

function rosen(x::Vector)
    (x[1])^2 + exp(x[2])
end

default(size=(600,600), fc=:heat)
x, y = -120.:0.5:120., -10.:0.1:10.
z = Surface((x,y)->rosen([x,y]), x, y)

surface(x,y,z)

#scatter3d!((0.,0.,rosen([0.,0.])))
x=-3:0.1:3
plot(x,x.^2, linewidth=4, linecolor=:orange)
plot!(x,x->3.0, linewidth=2)
