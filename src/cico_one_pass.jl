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
    loss_tol::Float64 = 1e-3, # we have no ideas how to implement loss_tol, currently it is in `inequality_constraint!`
    # better results for LN_NELDERMEAD, :LD_MMA, :LD_SLSQP, :LD_CCSAQ
    # worse results for :LN_COBYLA, :LN_PRAXIS
    # errors for :LN_BOBYQA, :LN_SBPLX, :LN_NEWUOA
    local_alg::Symbol = :LN_NELDERMEAD,
    # options for local fitter :max_iter
    max_iter::Int = 10^5,
    loss_grad::Union{Function, Symbol} = :EMPTY,
    kwargs...
    )
    # dim of the theta vector
    n_theta = length(theta_init)

    # checking arguments
    # when using :LN_NELDERMEAD initial parameters should not be zero
    if local_alg == :LN_NELDERMEAD
        zeroParameter = [ isapprox(theta_init[i], 0., atol=1e-2) for i in 1:n_theta]
        if any(zeroParameter)
            @warn "Close-to-zero parameters found when using :LN_NELDERMEAD."
            show(findall(zeroParameter))
        end
    end
    
    # checking loss_grad
    is_gradient = occursin(r"^LD_", String(local_alg))
    if loss_grad == :EMPTY && is_gradient
        throw(ArgumentError("`loss_grad` must be set for gradient local fitter `$(local_alg)`"))
    end

    # optimizer
    local_opt = Opt(local_alg, n_theta)
    ftol_abs!(local_opt, scan_tol) #ftol_abs

    # flags to analyze fitting stop
    out_of_bound::Bool = false
    #constraints_func_error = false

    function constraints_func(x, g) # testing grad methods
        # function constraints_func(x) # testing grad methods    
        # this part is necessary to understand the difference between
        # "stop out of bounds" and "stop because of function call error"
        loss = try
            loss_func(x)
        catch e
            @warn "Error when call loss_func($x)"
            throw(e)
        end
        
        if (loss < 0.) && (scan_func(x) > scan_bound)
            out_of_bound = true
            throw(NLopt.ForcedStop())
        elseif length(g) > 0
          if loss_grad == :AUTODIFF
            ForwardDiff.gradient!(g, loss_func, x)
          elseif loss_grad == :FINITE
            Calculus.finite_difference!(loss_func,x,g,:central)
          end
        end
        
        return loss
    end

    function obj_func(x, g)
        if length(g) > 0
          if loss_grad == :AUTODIFF
            ForwardDiff.gradient!(g, scan_func, x)
          elseif loss_grad == :FINITE
            Calculus.finite_difference!(scan_func, x, g, :central)
          end
        end
        return scan_func(x)
    end

    # constrain optimizer
    opt = Opt(:LN_AUGLAG, n_theta)
    ftol_abs!(opt, scan_tol)

    max_objective!(
        opt,
        obj_func
        )
        
    local_optimizer!(opt, local_opt)
    maxeval!(opt, max_iter)

    # inequality constraints
    inequality_constraint!(
        opt,
        constraints_func,
        loss_tol
    )

    # version 1: internal :LN_AUGLAG box constrains
    opt.lower_bounds = [tb[1] for tb in theta_bounds]
    opt.upper_bounds = [tb[2] for tb in theta_bounds]

    # version 2: creating lower and upper bounds as general constraints
    #= 
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
    =#

    # start optimization
    (optf, optx, ret) = optimize(opt, theta_init)

    if ret == :FORCED_STOP && !out_of_bound
        pp = ProfilePoint[]
        res = (nothing, pp, :LOSS_ERROR_STOP)
    elseif ret == :MAXEVAL_REACHED
        pp = ProfilePoint[]
        res = (nothing, pp, :MAX_ITER_STOP)
    elseif (ret == :FORCED_STOP || ret == :FAILURE) && out_of_bound # successfull result
        pp = ProfilePoint[]
        res = (nothing, pp, :SCAN_BOUND_REACHED)
    elseif ret == :FTOL_REACHED # successfull result
        loss = loss_func(optx)
        pp = [ ProfilePoint(optf, loss, optx, ret, nothing) ]
        res = (optf, pp, :BORDER_FOUND_BY_SCAN_TOL)
    else
        # this part is reached in case of :FAILURE return code
        # we have seen it in case of gradient methods
        pp = ProfilePoint[]
        res = (nothing, pp, :UNKNOWN_STOP)
    end

    return res
end # get_right_endpoint

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
    kwargs... # for get_right_endpoint
    )
    # checking arguments
    if theta_num > length(theta_init)
        throw(DomainError(theta_num, "theta_num exceed theta dimention"))
    end

    scan_func(theta::Vector) = theta[theta_num]

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
        kwargs...
    )
end

# struct OutOfBoundException <: Exception end
