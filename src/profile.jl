using NLopt

"""
Returns profile function for selected parameter component"
"""
function profile(
    theta_init::Vector{Float64},
    theta_num::Int,
    loss_func::Function;

    skip_optim::Bool = false,
    scale::Vector{Symbol} = fill(:direct, length(theta_init)),
    theta_bounds::Vector{Tuple{Float64,Float64}} = unscaling.(
        fill((-Inf, Inf), length(theta_init)),
        scale
        ),
    loss_tol::Float64 = 1e-3,
    local_alg::Symbol = :LN_NELDERMEAD,
    max_iter::Int = 10^5,
    kwargs... # options for local fitter
    )
    theta_length = length(theta_init)
    # set indexes
    indexes_rest = [i for i in 1:theta_length]
    deleteat!(indexes_rest, theta_num)
    # set bounds
    lb = [theta_bounds[i][1] for i in indexes_rest]
    ub = [theta_bounds[i][2] for i in indexes_rest]

    if skip_optim || theta_length == 1 # if profile == loss_func
        return (x::Float64; theta_init_i::Vector{Float64} = theta_init) -> begin
            theta_full = copy(theta_init)
            splice!(theta_full, theta_num, x)
            loss = loss_func(theta_full)
            # return
            ProfilePoint(
                loss,
                theta_full,
                :OPTIMIZATION_SKIPPED
            )
        end
    else
        out = (x::Float64; theta_init_i::Vector{Float64} = theta_init) -> begin
            # get init of rest component
            theta_init_rest = copy(theta_init_i)
            deleteat!(theta_init_rest, theta_num)
            # get rest of loss_func
            loss_func_rest = (theta_rest::Array{Float64, 1}, g::Array{Float64, 1}) -> begin
                theta_full = copy(theta_rest)
                splice!(theta_full, theta_num:(theta_num-1), x)
                loss_func(theta_full)
            end
            # set optimizer
            opt = Opt(local_alg, length(theta_init_rest))
            min_objective!(opt, loss_func_rest)
            ftol_abs!(opt, loss_tol)
            lower_bounds!(opt, lb)
            upper_bounds!(opt, ub)
            maxeval!(opt, max_iter)
            # start optimization
            (loss, theta_opt, ret) = optimize(opt, theta_init_rest)
            splice!(theta_opt, theta_num:(theta_num-1), x)
            # return
            ProfilePoint(
                loss,
                theta_opt,
                ret
            )
        end
        return out
    end
end
