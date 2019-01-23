using RecipesBase

"""
    using Plots
    plotly()
    plot(pi::ParamInterval)

Plots profile `L(theta)` for parameter `theta_num`,
`identifiability level`, `identifiability interval`.
Use `update_profile_points!(pi::ProfileInterval)` function to refine
profile points and make your plot more smooth
"""
@recipe function f(pi::ParamInterval)
    # get points for plot
    xs, ys = get_grid(pi)
    xlabel --> "theta"
    ylabel --> "L(theta)"

    if length(xs) > 1
        @series begin
            label --> "Theta parameter profile points"
            seriestype --> :line
            linestyle --> :dot
            markershape --> :circle
            markercolor --> :green
            line := (2.2,:green)
            (xs, ys)
        end
    end

    # critical level subplot
    @series begin
        label --> "Identifiability level"
        seriestype --> :hline
        line := (2.5, :purple)
        [pi.input.loss_crit]
    end

    if (pi.result[1].value != nothing) || (pi.result[2].value != nothing)
        @series begin
            label --> "Identifiability inteval"
            seriestype --> :vline
            line := (2.1, :pink)
            [pi.result[1].value != nothing ? pi.result[1].value : NaN,
            pi.result[2].value != nothing ? pi.result[2].value : NaN]

            #[pi.result[1].value,pi.result[2].value]
        end
    end

    @series begin
        label --> "Initial point"
        seriestype --> :scatter
        markershape --> :diamond
        markercolor --> :orange
        [(pi.input.theta_init[pi.input.theta_num], pi.loss_init)]
    end

end


function get_grid(pi::ParamInterval)
    xs = Vector{Float64}()
    ys = Vector{Float64}()
    #ns = Vector{String}()

    ep_start = pi.input.theta_init[pi.input.theta_num]
    push!(xs, ep_start)
    push!(ys, pi.loss_init)
    #push!(ns, "1")

    for ep in pi.result
        x_pps, y_pps = get_pps(ep)
        append!(xs, x_pps)
        append!(ys, y_pps)
        #append!(ns, n_pps)
    end
    return (xs, ys)
end

function get_pps(ep::EndPoint)
    l = length(ep.profilePoints)

    xs = Vector{Float64}(undef, l)
    ys = Vector{Float64}(undef, l)
    #ns = Vector{String}(undef, l)

    for i in 1:l
        xs[i] = ep.profilePoints[i].value
        ys[i] = ep.profilePoints[i].loss
        #ns[i] = ep.profilePoints[i].counter != nothing ?
        #    string(ep.profilePoints[i].counter) : ""
    end
    (xs,ys)
end

"""
    update_profile_points!(pi::ParamInterval)

Refines profile points to make your plot more smooth. Internally uses
`adapted_grid` to compute additional profile points.
See `PlotUtils.adapted_grid`
## Arguments
* `max_recursions`: how many times each interval is allowed to
be refined (default: 2).
"""
function update_profile_points!(pi::ParamInterval; max_recursions::Int = 2)
    for ep in pi.result
        if ep.status != :SCAN_BOUND_REACHED
            ep_start = pi.input.theta_init[pi.input.theta_num]
            update_profile_endpoint!(
                ep.profilePoints,
                ep.direction == :left ?
                    (ep.value, ep_start) :
                    (ep_start, ep.value),
                pi.input.theta_init,
                pi.input.theta_num,
                pi.input.loss_func,
                pi.input.local_alg,
                pi.input.theta_bounds,
                pi.input.loss_tol;
                max_recursions = max_recursions
            )
        else
            @info "$(string(ep.direction)) - a half-open interval"
        end
    end
    return nothing
end


function update_profile_endpoint!(
    pps_arr::Vector{ProfilePoint},
    interval::Tuple{Float64,Float64},
    init_params::Vector{Float64},
    id::Int64,
    loss_func::Function,
    fit_alg::Symbol,
    bounds::Vector{Tuple{Float64,Float64}},
    tol::Float64;
    max_recursions::Int = 2
)
    params = copy(init_params)

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

            loss = loss_func(minx)
            push!(pps_arr, ProfilePoint(x, loss, minx, ret, nothing))
            return loss
        end
    end

    # adapted_grid
    adapted_grid2(profile_func,interval; max_recursions = max_recursions)
    return nothing
end

# deprecated
function get_adapted_grid(
    interval::Tuple{Float64,Float64},
    init_params::Vector{Float64},
    id::Int64,
    loss_func::Function,
    fit_alg::Symbol,
    bounds::Vector{Tuple{Float64,Float64}},
    tol::Float64;
    max_recursions::Int = 1
)
    params = copy(init_params)
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
