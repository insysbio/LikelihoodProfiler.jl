using NLopt

"""
    function profile(
        theta_init::Vector{Float64},
        theta_num::Int,
        loss_func::Function;

        skip_optim::Bool = false,
        theta_bounds::Vector{Tuple{Float64,Float64}} = fill((-Inf, Inf), length(theta_init)),
        local_alg::Symbol = :LN_NELDERMEAD,
        ftol_abs::Float64 = 1e-3,
        maxeval::Int = 10^5,
        kwargs...
        )
Returns profile function for selected parameter component.
"""
function profile(
    theta_init::Vector{Float64},
    theta_num::Int,
    loss_func::Function;

    skip_optim::Bool = false,
    theta_bounds::Vector{Tuple{Float64,Float64}} = fill((-Inf, Inf), length(theta_init)),
    # fit alg args
    local_alg::Symbol = :LN_NELDERMEAD,
    ftol_abs::Float64 = 1e-3,
    kwargs...
    )
    theta_length = length(theta_init)

    # set indexes
    indexes_rest = [i for i in 1:theta_length]
    deleteat!(indexes_rest, theta_num)
    # set bounds
    lb = [theta_bounds[i][1] for i in indexes_rest]
    ub = [theta_bounds[i][2] for i in indexes_rest]

    if skip_optim || theta_length == 1 # if profile == loss_func
        return (x::Float64; theta_init_i::Vector{Float64} = theta_init, maxeval::Int = 10^5) -> begin
            theta_full = copy(theta_init)
            splice!(theta_full, theta_num, x)
            loss = loss_func(theta_full)
            # return
            ProfilePoint(
                x,
                loss,
                theta_full,
                :OPTIMIZATION_SKIPPED,
                1 # counter, only one call
            )
        end
    else
        # set optimizer
        opt = Opt(local_alg, theta_length - 1)
        ftol_abs!(opt, ftol_abs)
        lower_bounds!(opt, lb)
        upper_bounds!(opt, ub)
        # profile function
        return (x::Float64; theta_init_i::Vector{Float64} = theta_init, maxeval::Int = 10^5) -> begin
            # to count loss function calls inside profile
            counter::Int = 0
            # get init of rest component
            theta_init_rest = copy(theta_init_i)
            deleteat!(theta_init_rest, theta_num)
            # get rest of loss_func
            loss_func_rest = (theta_rest::Array{Float64, 1}, g::Array{Float64, 1}) -> begin
                theta_full = copy(theta_rest)
                splice!(theta_full, theta_num:(theta_num-1), x)
                counter += 1 # update counter
                loss_func(theta_full)
            end
            # set optimizer
            min_objective!(opt, loss_func_rest)
            maxeval!(opt, maxeval)
            # start optimization
            (loss, theta_opt, ret) = optimize(opt, theta_init_rest)
            splice!(theta_opt, theta_num:(theta_num-1), x)
            # return
            ProfilePoint(
                x,
                loss,
                theta_opt,
                ret,
                counter
            )
        end
    end
end
