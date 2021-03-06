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
    out_of_bound::Bool = false

    function constraints_func(x, g) # testing grad methods
    #function constraints_func(x) # testing grad methods    
        # this part is necessary to understand the difference between
        # "stop out of bounds" and "stop because of function call error"
        try
            loss = loss_func(x)
        catch e
            @warn "Error when call loss_func($x)"
            throw(e)
        end
        #println("constr")
        #@show (x,g)
        #@show (loss)
        #@show scan_bound
        if (loss < 0.) && (scan_func(x) > scan_bound)
            out_of_bound = true
            throw(NLopt.ForcedStop())
        #elseif isapprox(loss, 0., atol=loss_tol)
            #@warn "loss_tol reached... but..."
            #return loss
        elseif length(g) > 0
            try ForwardDiff.gradient!(g, loss_func, x)
            catch e
                @show e 
            end
            #Calculus.finite_difference!(loss_func,x,g,:central)
            #g .= Zygote.gradient(loss_func,x)[1]
        end
        
        return loss
    end

    function obj_func(x,g)
        #println("obj")
        #@show (x,g)
        if length(g) > 0
            try ForwardDiff.gradient!(g, scan_func, x)
            catch e
                @show e 
            end
            #Calculus.finite_difference!(scan_func,x,g,:central)
            #g .= Zygote.gradient(scan_func,x)[1]
        end
        
        scan_func(x)
    end

    # constrain optimizer
    opt = Opt(:LN_AUGLAG, n_theta)
    ftol_abs!(opt, scan_tol)

    max_objective!(
        opt,
        obj_func
        #NLoptAdapter(scan_func,theta_init)
        )
        
    local_optimizer!(opt, local_opt)
    maxeval!(opt, max_iter)

    # inequality constraints
    inequality_constraint!(
        opt,
        #NLoptAdapter(constraints_func,theta_init),
        constraints_func,
        loss_tol
    )
    opt.lower_bounds = [tb[1] for tb in theta_bounds]
    opt.upper_bounds = [tb[2] for tb in theta_bounds]
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
        loss = loss_func(optx)
        pp = [ ProfilePoint(optf, loss, optx, ret, nothing) ]
        res = (optf, pp, :BORDER_FOUND_BY_SCAN_TOL)
    else
        # this part is not normally reached, just for case
        #throw(ErrorException("No interpretation of the optimization results."))
        # do not throw
        pp = ProfilePoint[]
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
        kwargs... # options for local fitter
    )
end
