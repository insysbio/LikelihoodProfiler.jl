
# evaluate right bound of scan_val
function get_right_endpoint(
    theta_init::Vector{Float64}, # initial point of parameters
    scan_loss_func::Function, # function returns tuple of scan value, loss value: (h, lambda)
    method::Val{:CICO_ONE_PASS}; # function works only for method ONE_PASS;

    theta_bounds::Vector{Tuple{Float64,Float64}} = fill(
        (-Inf, Inf), length(theta_init)
        ),
    scan_bound::Float64 = 9.0,
    scan_tol::Float64 = 1e-3,
    loss_tol::Float64 = 1e-3, # i do not know how to implement it it
    # good results for :LN_NELDERMEAD, :LN_COBYLA, :LN_PRAXIS,
    # errors in :LN_BOBYQA, :LN_SBPLX, :LN_NEWUOA
    box_theta::Bool = false, # EXPERMENTAL, if true scan_loss_func cannot be caculated out of theta_bounds
    local_alg::Symbol = :LN_NELDERMEAD,
    # options for local fitter :max_iter
    max_iter::Int = 10^5,
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

    # Constraints function: loss_val <= 0
    out_of_bound::Bool = false
    function constraints_func(x::Vector{Float64}, g)
        try
            boxed_theta = box_theta ? boxing(x, theta_bounds) : x
            (scan_val, loss_val) = scan_loss_func(boxed_theta)          # call 1
        catch e
            @warn "Error when call scan_loss_func($x)"
            throw(e)
        end

        # this part is necessary to understand the difference between
        # "stop out of bounds" and "stop because of function call error"
        if (loss_val < 0.) && (scan_val > scan_bound)
            out_of_bound = true
            throw(ForcedStop("Out of the scan bound but in ll constraint."))
        #elseif isapprox(loss_val, 0., atol=loss_tol)
            #@warn "loss_tol reached... but..."
            #return loss_val
        end

        return loss_val
    end

    # condition for scan_val
    opt = Opt(:LN_AUGLAG, n_theta)
    ftol_abs!(opt, scan_tol)
    max_objective!(
        opt,
        function(x::Vector{Float64}, g)
            boxed_theta = box_theta ? boxing(x, theta_bounds) : x
            (scan_val, loss_val) = scan_loss_func(boxed_theta)          # call 2
            scan_val
        end
        )
    local_optimizer!(opt, local_opt)
    maxeval!(opt, max_iter)

    # inequality constraints
    inequality_constraint!(
        opt,
        constraints_func,
        loss_tol
    )
    [ inequality_constraint!(
        opt,
        (x, g) -> x[i] - theta_bounds[i][2],
        0.
    ) for i in 1:n_theta ]
    [ inequality_constraint!(
        opt,
        (x, g) -> theta_bounds[i][1] - x[i],
        0.
    ) for i in 1:n_theta ]

    # start optimization: (max scan_val, optimal params, code)
    (optf, optx, ret) = optimize(opt, theta_init)

    if (ret == :FORCED_STOP && !out_of_bound)
        pp = ProfilePoint[]
        res = (nothing, pp, :LOSS_ERROR_STOP)
    elseif ret == :MAXEVAL_REACHED
        pp = ProfilePoint[]
        res = (nothing, pp, :MAX_ITER_STOP)
    elseif (ret == :FORCED_STOP && out_of_bound) # successfull result
        pp = ProfilePoint[]
        res = (nothing, pp, :SCAN_BOUND_REACHED)
    elseif ret == :FTOL_REACHED # successfull result
        boxed_theta = box_theta ? boxing(optx, theta_bounds) : optx
        (scan_val, loss_val) = scan_loss_func(boxed_theta)              # call 3
        pp = [ ProfilePoint(optf, loss_val, boxed_theta, ret, nothing) ]
        res = (optf, pp, :BORDER_FOUND_BY_SCAN_TOL)
    else
        # this part is not normally reached, just for case
        throw(ErrorException("No interpretation of the optimization results: $ret"))
        # do not throw
        #pp = ProfilePoint[]
        #res = (nothing, pp, :UNKNOWN_STOP)
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
    box_theta::Bool = false, # EXPERMENTAL, if true loss_func cannot be caculated out of theta_bounds
    local_alg::Symbol = :LN_NELDERMEAD,
    kwargs... # options for local fitter
    )
    # checking arguments
    if theta_num > length(theta_init)
        throw(DomainError(theta_num, "theta_num exceed theta dimention"))
    end

    function scan_loss_func(theta::Vector{Float64})
        (theta[theta_num], loss_func(theta))
    end

    get_right_endpoint(
        theta_init,
        scan_loss_func,
        method;

        theta_bounds = theta_bounds,
        scan_bound = scan_bound,
        scan_tol = scan_tol,
        loss_tol = loss_tol,
        box_theta = box_theta,
        local_alg = local_alg,
        kwargs... # options for local fitter
    )
end

function boxing(
    theta::Vector{Float64},
    theta_bounds::Vector{Tuple{Float64,Float64}}
    )
    #checking arguments
    if length(theta) !== length(theta_bounds)
        throw(DimensionMismatch("x and box dimention mismatch: "*string(length(theta))*" vs "*string(length(theta_bounds))))
    end

    one_component = (i::Int64) ->
        if theta[i] < theta_bounds[i][1]
            theta_bounds[i][1]
        elseif theta[i] > theta_bounds[i][2]
            theta_bounds[i][2]
        else
            theta[i]
    end

    [one_component(i) for i in 1:length(theta)]
end
