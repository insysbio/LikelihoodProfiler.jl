using NLopt

# evaluate right bound of scan_func
function get_right_endpoint(
    theta_init::Vector{Float64}, # initial point of parameters
    scan_func::Function, # h(theta) function for predictions or parameters
    loss_func::Function, # lambda(theta) - labmbda_min - delta_lambda
    method::Val{:CICO_ONE_PASS}; # function works only for method ONE_PASS;

    theta_bounds::Vector{Tuple{Float64,Float64}} = fill(
        (-Inf, Inf), length(theta_init)
        ),
    scan_bound::Float64 = 9.0,
    scan_tol::Float64 = 1e-3,
    loss_tol::Float64 = 1e-3, # i do not know how to use it
    # good results in :LN_NELDERMEAD, :LN_COBYLA, :LN_PRAXIS,
    # errors in :LN_BOBYQA, :LN_SBPLX, :LN_NEWUOA
    local_alg::Symbol = :LN_NELDERMEAD,
    # options for local fitter :max_iter
    max_iter::Int = 10^5,
    ftol_abs::Float64 = 1e-3,
    kwargs...
    )
    # dim of the theta vector
    n_theta = length(theta_init)

    # checking arguments
    # methods which are not supported
    if local_alg in [:LN_BOBYQA, :LN_SBPLX, :LN_NEWUOA]
        @warn "Using local_alg = :"*String(local_alg)*" may result in wrong output."
    end
    # when using :LN_NELDERMEAD initial parameters should not be zero
    if local_alg == :LN_NELDERMEAD
        zeroParameter = [ isapprox(theta_init[i], 0., atol=1e-2) for i in 1:n_theta]
        if any(zeroParameter)
            @warn "Close-to-zero parameters found when using :LN_NELDERMEAD."
            show(findall(zeroParameter))
        end
    end

    # optimizer
    local_opt = Opt(local_alg, n_theta)
    ftol_abs!(local_opt, scan_tol) #ftol_abs

    # Constraints function
    function constraints_func(x, g)
        loss = loss_func(x)
        if (loss < 0.) && (scan_func(x) > scan_bound)
            throw(ForcedStop("Out of the scan bound but in ll constraint."))
        #elseif isapprox(loss, 0., atol=loss_tol)
        #    @warn "loss_tol reached... but..."
        #    return loss
        else
            return loss
        end
    end

    # constrain optimizer
    opt = Opt(:LN_AUGLAG, n_theta)
    ftol_abs!(opt, scan_tol)
    max_objective!(
        opt,
        (x, g) -> scan_func(x)
        )
    lb = [theta_bounds[i][1] for i in 1:n_theta] # minimum.(theta_bounds)
    ub = [theta_bounds[i][2] for i in 1:n_theta] # maximum.(theta_bounds)
    lower_bounds!(opt, lb)
    upper_bounds!(opt, ub)
    local_optimizer!(opt, local_opt)
    maxeval!(opt, max_iter)

    # inequality constraints
    inequality_constraint!(
        opt,
        constraints_func,
        loss_tol
    )

    # start optimization
    (optf, optx, ret) = optimize(opt, theta_init)

    if ret == :FORCED_STOP
        pp = []
        res = (nothing, pp, :SCAN_BOUND_REACHED)
    elseif ret == :FTOL_REACHED
        loss = loss_func(optx)
        pp = [ ProfilePoint(optf, loss, optx, ret, nothing) ]
        res = (optf, pp, :BORDER_FOUND_BY_SCAN_TOL)
    else
        pp = []
        res = (nothing, pp, :UNKNOWN_STOP)
    end

    return res
end # of bound_right

function get_right_endpoint(
    theta_init::Vector{Float64}, # initial point of parameters
    theta_num::Int, # number of parameter to scan
    loss_func::Function, # lambda(theta) - labmbda_min - delta_lambda
    method::Val{:CICO_ONE_PASS}; # function works only for method ONE_PASS;

    theta_bounds::Vector{Tuple{Float64,Float64}} = fill(
        (-Inf, Inf), length(theta_init)
        ),
    scan_bound::Float64 = 9.0,
    scan_tol::Float64 = 1e-3,
    loss_tol::Float64 = 1e-3,
    local_alg::Symbol = :LN_NELDERMEAD,
    kwargs... # options for local fitter
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
        scan_bound = scan_bound,
        scan_tol = scan_tol,
        loss_tol = loss_tol,
        local_alg = local_alg,
        kwargs... # options for local fitter
    )
end
