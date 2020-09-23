"""
Implementation of an augmented Lagrangian method based on NLopt package.
"""

function auglag(
    ::Val{NLopt},
    theta_init::Vector{Float64}, # initial point of parameters
    scan_func, # function returns scan value
    loss_func; # function returns loss value
    local_alg::LN_NELDERMEAD
)
    # optimizer

    local_opt = Opt(local_alg, n_theta)
    ftol_abs!(local_opt, scan_tol) #ftol_abs

    # Constraints function: loss_val <= 0
    out_of_bound = false
    function constraints_func(x, g)
        if length(g)>0
            Calculus.finite_difference!(loss_func,x,g,:central)
            #ForwardDiff.gradient!(g, loss_func, x)
        end
        try
            scan_val = scan_func(x)
            loss_val = loss_func(x)
        catch e
            msg = e.msg
            @warn "Error when call loss_func($x) for loss_val. $msg"
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
        function(x, g)

            if length(g)>0
                Calculus.finite_difference!(scan_func,x,g,:central)
                #ForwardDiff.gradient!(g,scan_func, x)
            end
            scan_func(x)
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

    function left_bound_func(x,grad,theta_bounds,i)
        if length(grad)>0
            #grad .= zeros(length(grad))
            grad[i] = -1.0
        end
        theta_bounds[i][1] - x[i]
    end

    function right_bound_func(x,grad,theta_bounds,i)
        if length(grad)>0
            #grad .= zeros(length(grad))
            grad[i] = 1.0
        end
        x[i] - theta_bounds[i][2]
    end

    [ inequality_constraint!(
        opt,
        (x, g) -> right_bound_func(x,g,theta_bounds,i),
        0.
    ) for i in 1:n_theta ]

    [ inequality_constraint!(
        opt,
        (x, g) -> left_bound_func(x,g,theta_bounds,i),
        0.
    ) for i in 1:n_theta ]

    # start optimization: (max scan_val, optimal params, code)

    (optf, optx, ret) = optimize(opt, theta_init)
end


"""
    auglag(nlp)

Implementation of an augmented Lagrangian method. The following keyword parameters can be passed:
- μ: Starting value of the penalty parameter (default: 10.0)
- atol: Absolute tolerance used in dual feasibility measure (default: 1e-8)
- rtol: Relative tolerance used in dual feasibility measure (default: 1e-8)
- ctol: (Absolute) tolerance used in primal feasibility measure (default: 1e-8)
- max_iter: Maximum number of iterations (default: 1000)
- max_time: Maximum elapsed time in seconds (default: 30.0)
- max_eval: Maximum number of objective function evaluations (default: 100000)
- subsolver_logger: Logger passed to `tron` (default: NullLogger)
"""

function auglag(
    x0::Vector{T}, # initial point of parameters
    scan_func, # function returns scan value
    loss_func; # function returns loss value
    subsolver::Symbol = :LN_NELDERMEAD,
    μ::Real = eltype(x0)(10.0),
    max_iter::Int = 1000,
    max_time::Real = 30.0,
    max_eval::Int=100000,
    atol::Real = 1e-8,
    rtol::Real = 1e-8,
    ctol::Real = 1e-8,
    stopval::Real = eltype(x0)(1e9)
) where T
    T = eltype(x0)


end
