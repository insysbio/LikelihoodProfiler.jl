using Plots

"""
    using Plots
    plot(pi::ParamInterval)

Plots profile for parameter `theta_num`.
Computes `adapted_grid` for :CICO_ONE_PASS method
See also: `Plots.adapted_grid`

# Arguments
"""
@recipe function f(pi::ParamInterval)
    xs = Vector{Float64}()
    ys = Vector{Float64}()
    for ep in pi.result
        x_pps, y_pps = get_pps(ep)
        append!(xs, x_pps)
        append!(ys, y_pps)

        if pi.method == :CICO_ONE_PASS && ep.status == :BORDER_FOUND_BY_SCAN_TOL
            ep_start = pi.input.theta_init[pi.input.theta_num]
            x_ag, y_ag = get_adapted_grid(
                ep.direction == :left ?
                    (ep.value, ep_start) :
                    (ep_start, ep.value),
                pi.input.theta_init,
                pi.input.theta_num,
                pi.input.loss_func,
                pi.input.local_alg,
                pi.input.theta_bounds,
                pi.input.loss_tol
            )
            append!(xs, x_ag)
            append!(ys, y_ag)
        end

    end
    @series begin
        seriestype --> :line
        markershape --> :circle
        line := (2.2,:green)
        (xs, ys)
    end

    @series begin
        subplot := 1
        line := (:2.5, :red)
        seriestype --> :hline
        [pi.input.loss_crit]
    end
    #seriestype  :=  :path
end

function get_pps(ep::EndPoint)
    pps = ep.profilePoints
    xs = [pps[i].value for i in 1:length(pps)]
    ys = [pps[i].loss for i in 1:length(pps)]
    (xs,ys)
end

function get_adapted_grid(
    interval::Tuple{Float64,Float64},
    params::Vector{Float64},
    id::Int64,
    loss_func::Function,
    fit_alg::Symbol,
    bounds::Vector{Tuple{Float64,Float64}},
    tol::Float64;
    max_recursions::Int = 2
)
    # bounds
    println(interval)
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
