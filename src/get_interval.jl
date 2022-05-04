"""
    struct ParamIntervalInput
        theta_init::Vector{Float64} # initial parameters vector
        theta_num::Int # number of the parameter for identifiability analysis
        scan_func::Function # scan function
        loss_func::Function # loss function
        loss_crit::Float64 # loss function maximum value, "identifiability level"
        scale::Vector{Symbol}
        theta_bounds::Vector{Tuple{Float64, Float64}} # search bounds for id parameter
        scan_bounds::Tuple{Float64,Float64}
        scan_tol::Float64 # fitting tolerance for local optimizer (default - `1e-3`)
        loss_tol::Float64 # constraints tolerance
        local_alg::Symbol # local fitting algorithm (default - `:LN_NELDERMEAD`)
        fitter_options::Any
    end

Structure storing input data for parameter identification
"""
struct ParamIntervalInput
    theta_init::Vector{Float64} # initial parameters vector
    theta_num::Int # number of the parameter for analysis
    scan_func::Function
    loss_func::Function # loss function
    loss_crit::Float64 # loss function maximum value, "identifiability level"
    scale::Vector{Symbol}
    theta_bounds::Vector{Tuple{Float64, Float64}} # search bounds for id parameter
    scan_bounds::Tuple{Float64,Float64}
    scan_tol::Float64 # fitting tolerance for local optimizer (default - 1e-3)
    loss_tol::Float64 # constraints tolerance
    local_alg::Symbol # local fitting algorithm (default - :LN_NELDERMEAD)
    fitter_options::Any
end

"""
    struct ParamInterval
        input::ParamIntervalInput
        loss_init::Float64
        method::Symbol
        result::Tuple{EndPoint, EndPoint}
    end

Structure storing result of parameter identification
"""
struct ParamInterval
    input::ParamIntervalInput
    loss_init::Float64
    method::Symbol
    result::Tuple{EndPoint, EndPoint}
end

"""
    get_interval(
        theta_init::Vector{Float64},
        theta_num::Int,
        loss_func::Function,
        method::Symbol;

        loss_crit::Float64 = 0.0,
        scale::Vector{Symbol} = fill(:direct, length(theta_init)),
        theta_bounds::Vector{Tuple{Float64,Float64}} = unscaling.(
            fill((-Inf, Inf), length(theta_init)),
            scale
            ),
        scan_bounds::Tuple{Float64,Float64} = unscaling.(
            (-9.0, 9.0),
            scale[theta_num]
            ),
        scan_tol::Float64 = 1e-3,
        loss_tol::Float64 = 1e-3,
        local_alg::Symbol = :LN_NELDERMEAD,
        autodiff::Bool = true,
        kwargs...
        )
Computes confidence interval for single component `theta_num` of parameter vector.

## Return
`ParamInterval` structure storing all input data and estimated confidence interval.

## Arguments
- `theta_init`: starting values of parameter vector ``\\theta``. The starting values should not necessary be the optimum values of `loss_func` but `loss_func(theta_init)` should be lower than `loss_crit`.
- `theta_num`: index of vector component for identification: `theta_init(theta_num)`.
- `loss_func`: loss function ``\\Lambda\\left(\\theta\\right)`` for profile likelihood-based (PL) identification. Usually we use log-likelihood for PL analysis: ``\\Lambda( \\theta ) = - 2 ln\\left( L(\\theta) \\right)``.
- `method`: computational method to estimate confidence interval's endpoint. Currently the following methods are implemented: `:CICO_ONE_PASS`, `:LIN_EXTRAPOL`, `:QUADR_EXTRAPOL`.

## Keyword arguments
- `loss_crit`: critical level of loss function. Confidence interval's endpoint value is the intersection point of profile likelihood and `loss_crit` level.
- `scale`: vector of scale transformations for each parameters' component. Possible values: `:direct` (`:lin`), `:log`, `:logit`. This option can speed up the optimization, especially for wide `theta_bounds`. The default value is `:direct` (no transformation) for all parameters.
- `theta_bounds`: vector of tuple `(lower_bound, upper_bound)` for each parameter. Bounds define the ranges for possible parameter values. Default bounds are `(-Inf,Inf)`.
- `scan_bounds`: scan bounds tuple for `theta_num` parameter. Should be within the `theta_bounds` for `theta_num` parameter. Default is `(1e-9,1e9)`.
- `scan_tol`: Absolute tolerance for `theta_num` parameter used as termination criterion.  
- `loss_tol`: Absolute tolerance controlling `loss_func` closenes to `loss_crit` (termination criterion). Currently doesn't work for `:CICO_ONE_PASS` method because of limitation in `LN_AUGLAG` interface.
- `local_alg`: algorithm of optimization. Derivative-free and gradient-based algorithms form NLopt package. 
- `autodiff` : whether to use automatic differentiation with gradient-based algorithms. Default is `true`.
- `kwargs...`: the additional `method` specific keyword arguments.

"""
function get_interval(
    theta_init::Vector{Float64},
    theta_num::Int,
    loss_func::Function,
    method::Symbol;

    loss_crit::Float64 = 0.0,
    scale::Vector{Symbol} = fill(:direct, length(theta_init)),
    theta_bounds::Vector{Tuple{Float64,Float64}} = unscaling.(
        fill((-Inf, Inf), length(theta_init)),
        scale
        ),
    scan_bounds::Tuple{Float64,Float64} = unscaling.(
        (-9.0, 9.0),
        scale[theta_num]
        ),
    scan_tol::Float64 = 1e-3,
    loss_tol::Float64 = 1e-3,
    local_alg::Symbol = :LN_NELDERMEAD,
    kwargs... # other options for get_right_endpoint
    )
    # both endpoints
    endpoints = [get_endpoint(
        theta_init,
        theta_num,
        loss_func,
        method,
        [:left,:right][i]; # method
        loss_crit = loss_crit,
        scale = scale,
        theta_bounds = theta_bounds,
        scan_bound = scan_bounds[i],
        scan_tol = scan_tol,
        loss_tol = loss_tol,
        local_alg = local_alg,
        kwargs... # options for local fitter
        ) for i in 1:2]

    input = ParamIntervalInput(
        theta_init,
        theta_num,
        (x)->x[theta_num],
        loss_func,
        loss_crit,
        scale,
        theta_bounds,
        scan_bounds,
        scan_tol,
        loss_tol,
        local_alg,
        kwargs
        )

    ParamInterval(
        input,
        loss_func(theta_init),
        method,
        Tuple(endpoints)
        )
end


"""
    get_interval(
        theta_init::Vector{Float64},
        scan_func::Function,
        loss_func::Function,
        method::Symbol;

        loss_crit::Float64 = 0.0,
        scale::Vector{Symbol} = fill(:direct, length(theta_init)),
        theta_bounds::Vector{Tuple{Float64,Float64}} = unscaling.(
            fill((-Inf, Inf), length(theta_init)),
            scale
            ),
        scan_bounds::Tuple{Float64,Float64} = unscaling.(
            (-9.0, 9.0),
            :direct
            ),
        scan_tol::Float64 = 1e-3,
        loss_tol::Float64 = 1e-3,
        local_alg::Symbol = :LN_NELDERMEAD,
        autodiff::Bool = true,
        kwargs...
        )
Computes confidence interval for function of parameters `scan_func`.

## Return
`ParamInterval` structure storing all input data and estimated confidence interval.

## Arguments
- `theta_init`: starting values of parameter vector ``\\theta``. The starting values should not necessary be the optimum values of `loss_func` but `loss_func(theta_init)` should be lower than `loss_crit`.
- `scan_func`: scan function of parameters.
- `loss_func`: loss function ``\\Lambda\\left(\\theta\\right)`` the profile of which is analyzed. Usually we use log-likelihood for profile analysis in form ``\\Lambda( \\theta ) = - 2 ln\\left( L(\\theta) \\right)``.
- `method`: computational method to estimate confidence interval's endpoint. Currently supports only `:CICO_ONE_PASS` method.

## Keyword arguments
- `loss_crit`: critical level of loss function. Confidence interval's endpoint value is the intersection point of profile likelihood and `loss_crit` level.
- `scale`: vector of scale transformations for each parameters' component. Possible values: `:direct` (`:lin`), `:log`, `:logit`. This option can speed up the optimization, especially for wide `theta_bounds`. The default value is `:direct` (no transformation) for all parameters.
- `theta_bounds`: vector of tuple `(lower_bound, upper_bound)` for each parameter. Bounds define the ranges for possible parameter values. Default bounds are `(-Inf,Inf)`.
- `scan_bounds`: scan bounds tuple for `scan_func` values. Default is `(1e-9, 1e9)` .
- `scan_tol`: Absolute tolerance for `theta_num` parameter used as termination criterion.  
- `loss_tol`: Absolute tolerance controlling `loss_func` closenes to `loss_crit` (termination criterion). Currently doesn't work for `:CICO_ONE_PASS` method because of limitation in `LN_AUGLAG` interface.
- `local_alg`: algorithm of optimization. Derivative-free and gradient-based algorithms form NLopt package. 
- `autodiff` : whether to use automatic differentiation with gradient-based algorithms. Default is `true`.
- `kwargs...`: the additional `method` specific keyword arguments.
"""
function get_interval(
    theta_init::Vector{Float64},
    scan_func::Function,
    loss_func::Function,
    method::Symbol;

    loss_crit::Float64 = 0.0,
    scale::Vector{Symbol} = fill(:direct, length(theta_init)),
    theta_bounds::Vector{Tuple{Float64,Float64}} = unscaling.(
        fill((-Inf, Inf), length(theta_init)),
        scale
        ),
    scan_bounds::Tuple{Float64,Float64} = (-1e9, 1e9), # log scan bound is not implemented
    scan_tol::Float64 = 1e-3,
    loss_tol::Float64 = 1e-3,
    local_alg::Symbol = :LN_NELDERMEAD,
    kwargs... # other options for get_right_endpoint
    )
    # both endpoints
    endpoints = [get_endpoint(
        theta_init,
        scan_func,
        loss_func,
        method,
        [:left,:right][i]; # method
        loss_crit = loss_crit,
        scale = scale,
        theta_bounds = theta_bounds,
        scan_bound = scan_bounds[i],
        scan_tol = scan_tol,
        loss_tol = loss_tol,
        local_alg = local_alg,
        kwargs...
        ) for i in 1:2]

    input = ParamIntervalInput(
        theta_init,
        0,
        scan_func,
        loss_func,
        loss_crit,
        scale,
        theta_bounds,
        scan_bounds,
        scan_tol,
        loss_tol,
        local_alg,
        kwargs
        )

    ParamInterval(
        input,
        loss_func(theta_init),
        method,
        Tuple(endpoints)
        )
end