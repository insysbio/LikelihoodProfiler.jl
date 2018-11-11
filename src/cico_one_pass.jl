using NLopt

# evaluate right bound of scan_func
function get_right_endpoint(
    theta_init::Vector{Float64}, # initial point of parameters
    scan_func::Function, # h(theta) function for predictions or parameters
    loss_func::Function, # lambda(theta) - labmbda_min - delta_lambda
    method::Val{:CICO_ONE_PASS}; # function works only for method ONE_PASS;

    theta_bounds::Vector{Vector{Float64}} = fill(
        [-Inf, Inf], length(theta_init)
    ),
    scan_bound::Float64 = 9.0,
    local_alg::Symbol = :LN_NELDERMEAD,
    scan_tol::Float64 = 1e-3,
    # params_tol::Float64 = 1e-3, # not sure this required
    ll_tol::Float64 = 1e-3, # i do not know how to use it
    max_iter::Int64 = 10^5,
    kwargs...
)
    # dim of the theta vector
    n_theta = length(theta_init)

    # optimizer
    local_opt = Opt(local_alg, n_theta)
    ftol_abs!(local_opt, scan_tol)

    # Constraints function
    function constraints_func(x, g)
        loss = loss_func(x)
        if (loss < 0.) && (scan_func(x) > scan_bound)
            throw(ForcedStop("Out of the scan bound but in ll constraint."))
        else
            return loss
        end
    end

    # constrain optimizer
    opt = Opt(:LN_AUGLAG, n_theta)
    max_objective!(
        opt,
        (x, g) -> scan_func(x)
    )
    lb = minimum.(theta_bounds)
    ub = maximum.(theta_bounds)
    lower_bounds!(opt, lb)
    upper_bounds!(opt, ub)
    local_optimizer!(opt, local_opt)
    maxeval!(opt, max_iter)
    inequality_constraint!(
        opt,
        constraints_func,
        ll_tol
    )

    (optf, optx, ret) = optimize(opt, theta_init)

    if ret == :FORCED_STOP
        pp = []
        res = (scan_bound, pp, :SCAN_BOUND_REACHED)
    elseif ret == :FTOL_REACHED
        loss = loss_func(optx)
        pp = [ ProfilePoint(loss, optx, ret) ]
        res = (optf, pp, :BORDER_FOUND_BY_FTOL)
    else
        pp = []
        res = (scan_bound, pp, :UNKNOWN_STOP)
    end

    return res
end # of bound_right

# evaluate right or left endpoint parameter
function get_right_endpoint(
    theta_init::Vector{Float64}, # initial point of parameters
    theta_num::Int64, # number of parameter to scan
    loss_func::Function, # lambda(theta) - labmbda_min - delta_lambda
    method::Val{:CICO_ONE_PASS}; # function works only for method ONE_PASS;

    theta_bounds::Vector{Vector{Float64}} = fill(
        [-Inf, Inf], length(theta_init)
    ),
    kwargs...
)
    # checking arguments
    if theta_num > length(theta_init)
        throw(DomainError(theta_num, "theta_num exceed theta dimention"))
    end

    function scan_func(theta::Vector{Float64})
        theta[theta_num]
    end

    get_right_endpoint(
        theta_init,
        scan_func,
        loss_func,
        method;

        theta_bounds = theta_bounds,
        kwargs...
    )
end
