#using Plots

"""
    params_plot(params::Vector{Float64}, id::Int, loss_func::Function,
    interval::Tuple{Float64,Float64}; <keyword arguments>)

Computes `adapted_grid` for `loss_func` and `id` parameter values from the `interval`.
See also: `Plots.adapted_grid`

# Arguments
- `fit_alg::Symbol`: fitting algorithm (default `:LN_NELDERMEAD`).
- `bounds::Vector{Vector{Float64}}`: bound constraints for all parameters (default `[-Inf,Inf]`).
- `tol::Float64`: fitting tolerance (default `ftol_abs = 1e-3`).
- `max_recursions::Int`: how many times each interval is allowed to be refined (default `2`).
"""
function params_plot(
    params::Vector{Float64},
    id::Int,
    loss_func::Function,
    interval::Tuple{Float64,Float64};
    fit_alg::Symbol = :LN_NELDERMEAD,
    bounds::Vector{Vector{Float64}} = fill(
        [-Inf, Inf],
        length(params)
    ),
    tol::Float64 = 1e-3,
    max_recursions::Int = 2
)

    # bounds
    lb = minimum.(bounds)
    ub = maximum.(bounds)

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
    hline!([loss_crit], label="crit_level $loss_crit")
    scatter!(grid[1], grid[2], label="plot_grid")

    # return plot
    return gr
    =#

    # adapted_grid
    adapted_grid2(profile_func,interval; max_recursions = max_recursions)
end

function params_plot(
    params_interval::ParamInterval;
    delta::Float64 = 0.2, # magic number
    max_recursions::Int = 2
)
    params = params_interval.input.init_params
    id = params_interval.input.id
    loss_func = params_interval.input.loss_func
    interval = (
        params_interval.ret_codes[1] != :BOUNDS_REACHED ? params_interval.interval[1]-delta : 0.0,
        params_interval.ret_codes[2] != :BOUNDS_REACHED ? params_interval.interval[2]+delta : params_interval.input.scan_bound[2]
    )
    fit_alg = params_interval.input.local_alg
    bounds = params_interval.input.bounds
    tol = params_interval.input.losstol


    # bounds
    lb = minimum.(bounds)
    ub = maximum.(bounds)

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
    adapted_grid2(profile_func,interval; max_recursions = max_recursions)
end
