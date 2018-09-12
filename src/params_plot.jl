using Plots

"""
# Input:
        params - initial parameters vector
        id - id of the parameter for analysis
        maxf - loss function maximum value, "identifiability level"
        loss_func - loss function
        interval - interval for plot
        fit_alg - fitting algorithm (default - :LN_NELDERMEAD)
        bounds_params - ???
        tol - fitting tolerance (default - 1e-3)
        max_recursions - ???
# Return:
        parameter profile grid
"""
function params_plot(
    params::Vector{Float64},
    id::Int64,
    maxf::Float64,
    loss_func::Function,
    interval::Tuple{Float64,Float64};
    fit_alg::Symbol = :LN_NELDERMEAD,
    bounds_params::Vector{Vector{Float64}} = fill(
        [-Inf, Inf],
        length(params)
    ),
    tol::Float64 = 1e-3,
    max_recursions::Int64 = 2
)
    # cheking arguments
    (loss_func(params) > maxf) && throw(ArgumentError("Check params and maxf: loss_func(params) should be <= maxf"))

    # bounds
    lb = minimum.(bounds_params)
    ub = maximum.(bounds_params)

    # function to be used in adapted_grid
    function profile_func(x::Float64)
        if length(params) == 1
            return loss_func([x])
        else
            # ! params_intervals_one_side can be used here
            fit_params_func = (p,g) -> loss_func(p)
            opt = Opt(fit_alg, length(params))

            params[id] = x
            min_objective!(opt, fit_params_func)

            # excluding params[id] from optimization
            lb[id] = x
            ub[id] = x
            lower_bounds!(opt,lb)
            upper_bounds!(opt,ub)
            ftol_abs!(opt,tol)

            (loss,minx,ret) = optimize(opt,params)

            return loss
        end
    end
    # adapted_grid
    #=
    grid = adapted_grid(profile_func,interval; max_recursions = max_recursions)
    println("grid = $grid")
    gr = plot(grid[1], grid[2], label="param_$id");
    hline!([maxf], label="crit_level $maxf")
    scatter!(grid[1], grid[2], label="plot_grid")

    # return plot
    return gr
    =#

    # adapted_grid
    adapted_grid(profile_func,interval; max_recursions = max_recursions)
end
