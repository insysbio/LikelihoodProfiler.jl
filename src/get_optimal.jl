"""
        function get_optimal(
            theta_init::Vector{Float64}, # initial point of parameters
            loss_func::Function; # lambda(theta)

            scale::Vector{Symbol} = fill(:direct, length(theta_init)),
            theta_bounds::Vector{Tuple{Float64,Float64}} = unscaling.(
                fill((-Inf, Inf), length(theta_init)),
                scale
                ),
            scan_tol::Float64 = 1e-3,
            loss_tol::Float64 = 1e-3,
            local_alg::Symbol = :LN_NELDERMEAD,
            silent::Bool = false,
            max_iter::Int = 10^5
        )

    Provides the optimization routine using the interface similar to [`get_endpoint`](@ref).
    Currently it uses standard NLopt optimization but allows to use parameter scaling and autodiff method.

    ## Return
    [`ProfilePoint`](@ref) object where `value` and `loss` properties are equal to the optimal (minimal) `loss_func` value.

    ## Arguments

    - `theta_init`: starting values of parameter vector ``\\theta``.
    - `loss_func`: loss function ``\\Lambda\\left(\\theta\\right)`` for profile likelihood-based (PL) identification. Usually we use log-likelihood for PL analysis: ``\\Lambda( \\theta ) = - 2 ln\\left( L(\\theta) \\right)``.
    
    ## Keyword arguments
    - `scale`: vector of scale transformations for each parameters' component. Possible values: `:direct` (`:lin`), `:log`, `:logit`. This option can speed up the optimization, especially for wide `theta_bounds`. The default value is `:direct` (no transformation) for all parameters.
    - `theta_bounds`: vector of tuple `(lower_bound, upper_bound)` for each parameter. Bounds define the ranges for possible parameter values. Default bounds are `(-Inf,Inf)`.
    - `scan_tol`: Absolute tolerance for all component of theta used as termination criterion.  
    - `loss_tol`: Absolute tolerance controlling `loss_func`.
    - `local_alg`: algorithm of optimization. Derivative-free and gradient-based algorithms form NLopt package.
    - `silent` : Boolean argument declaring whether we display the optimization progress. Default is `false`
    - `autodiff` (experimental) : whether to use automatic differentiation with gradient-based algorithms. Default is `true`.
"""
function get_optimal(
    theta_init::Vector{Float64}, # initial point of parameters
    loss_func::Function; # lambda(theta)

    scale::Vector{Symbol} = fill(:direct, length(theta_init)),
    theta_bounds::Vector{Tuple{Float64,Float64}} = unscaling.(
        fill((-Inf, Inf), length(theta_init)),
        scale
        ),
    scan_tol::Float64 = 1e-3,
    loss_tol::Float64 = 1e-3,
    local_alg::Symbol = :LN_NELDERMEAD,
    silent::Bool = false,
    max_iter::Int = 10^5,
    loss_grad::Union{Function, Symbol} = :EMPTY
    # kwargs... # options for local fitter
)
    # dim of the theta vector
    n_theta = length(theta_init)

    # checking arguments
    # theta_bound[1] < theta_init < theta_bound[2]
    theta_init_outside_theta_bounds = .! [theta_bounds[i][1] < theta_init[i] < theta_bounds[i][2] for i in 1:n_theta]
    if any(theta_init_outside_theta_bounds)
        throw(ArgumentError("theta_init is outside theta_bound: $(findall(theta_init_outside_theta_bounds))"))
    end
    # 0 <= theta_bound[1] for :log
    less_than_zero_theta_bounds = (scale .== :log) .& [theta_bounds[i][1] < 0 for i in 1:n_theta]
    if any(less_than_zero_theta_bounds)
        throw(ArgumentError(":log scaled theta_bound min is negative: $(findall(less_than_zero_theta_bounds))"))
    end
    # 0 <= theta_bounds <= 1 for :logit
    less_than_zero_theta_bounds = (scale .== :logit) .& [theta_bounds[i][1] < 0 || theta_bounds[i][2] > 1 for i in 1:n_theta]
    if any(less_than_zero_theta_bounds)
        throw(ArgumentError(":logit scaled theta_bound min is outside range [0,1]: $(findall(less_than_zero_theta_bounds))"))
    end
    
    # checking arguments
    # when using :LN_NELDERMEAD initial parameters should not be zero
    if local_alg == :LN_NELDERMEAD
        zeroParameter = [ isapprox(theta_init[i], 0., atol=1e-2) for i in 1:n_theta]
        if any(zeroParameter)
            @warn "Close-to-zero parameters found when using :LN_NELDERMEAD."
            show(findall(zeroParameter))
        end
    end

    opt = Opt(local_alg, n_theta)
    opt.ftol_abs = loss_tol
    opt.xtol_abs = scan_tol

    # progress info
    prog = ProgressUnknown("Fitter counter:"; spinner=false, enabled=!silent, showspeed=true)
    count = 0
    supreme = nothing

    loss_func_prog = function(theta)
        loss_res = loss_func(theta)

        count += 1
        if typeof(supreme)==Nothing || loss_res < supreme
            supreme = loss_res
        end
        ProgressMeter.update!(prog, count, spinner="⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"; showvalues = [(:supreme,round(supreme; sigdigits=4))])
        
        return loss_res
    end

    # transforming loss
    theta_init_g = scaling.(theta_init, scale)
    loss_func_g = if loss_grad == :EMPTY
        function(theta_g, grad) # no gradient, only for gradient-free methods
            theta = unscaling.(theta_g, scale)
            
            return loss_func_prog(theta)
        end
    elseif isa(loss_grad, Function)
        function(theta_g, grad) # gradient is directly declared
            theta = unscaling.(theta_g, scale)

            # rescaling gradient function
            if length(grad) > 0
                loss_grad_theta = loss_grad(theta)
                for i in 1:n_theta
                    if scale[i] == :log
                        grad[i] = loss_grad_theta[i] * theta[i] * log(10.)
                    elseif scale[i] == :logit
                        grad[i] = loss_grad_theta[i] * theta[i] * (1. - theta[i]) * log(10.)
                    else
                        grad[i] = loss_grad_theta[i]
                    end
                end
            end

            return loss_func_prog(theta)
        end
    elseif loss_grad == :AUTODIFF
        function(theta_g, grad) # gradient is directly declared
            theta = unscaling.(theta_g, scale)
            ForwardDiff.gradient!(grad, loss_func_prog, theta_g)

            # rescaling gradient function
            loss_grad_theta = loss_grad(theta)
            for i in 1:n_theta
                if scale[i] == :log
                    grad[i] = loss_grad_theta[i] * theta[i] * log(10.)
                elseif scale[i] == :logit
                    grad[i] = loss_grad_theta[i] * theta[i] * (1. - theta[i]) * log(10.)
                else
                    grad[i] = loss_grad_theta[i]
                end
            end

            return loss_func_prog(theta)
        end
    end

    opt.min_objective = loss_func_g
    opt.maxeval = max_iter

    # version 1: internal :LN_AUGLAG box constrains
    theta_bounds_g = scaling.(theta_bounds, scale)
    opt.lower_bounds = [tb[1] for tb in theta_bounds_g]
    opt.upper_bounds = [tb[2] for tb in theta_bounds_g]

    (optf, optx_g, ret) = optimize(opt, theta_init_g)
    optx = unscaling.(optx_g, scale)

    pp = ProfilePoint(optf, loss_func(optx), optx, ret, count)

    return pp
end
