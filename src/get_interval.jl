"Structure storing input data for parameter interval calculation"
struct ParamIntervalInput
    theta_init::Vector{Float64} # initial parameters vector
    theta_num::Int # number of the parameter for analysis
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

"Structure storing result of parameter interval calculation"
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
        kwargs...
        )
Computes confidence interval for single component `theta_num` of parameter vector
and `loss_func` according to `loss_crit` level.

## Return
`ParamInterval` structure storing all input data and estimated confidence interval.

## Arguments
- `theta_init`: starting values of parameter vector ``\\theta``. The starting values is not necessary to be the optimum values for `loss_func` but it the value of `loss_func` must be lower than `loss_crit`.
- `theta_num`: number ``n`` of vector component to compute confidence interval ``\\theta^n``.
- `loss_func`: loss function ``\\Lambda\\left(\\theta\\right)`` the profile of which is analyzed. Usually we use log-likelihood for profile analysis in form ``\\Lambda( \\theta ) = - 2 ln\\left( L(\\theta) \\right)``
- `method`: computational method to evaluate interval endpoint. Currently the following methods are implemented: `:CICO_ONE_PASS`.

## Keyword arguments
- `loss_crit`: critical level of loss function.
- `scale`: vector of scale transformations for each component. Possible values: `:direct, :log, :logit`. This option can make optimization much more faster, especially for wide `theta_bounds`.
- `theta_bounds`: vector of bounds for each component in format `(left_border, right_border)`. The values outside the bound will be ignored.
- `scan_bounds`: vector of scan bound for `theta_num` component. It must be within the `theta_bounds` for the scanned component.
- `scan_tol`: Abolute tolerance of scanned component (stop criterion).
- `loss_tol`: *experimental*. Required tolerance of `loss_func` at `loss_crit`.
- `local_alg`: algorithm of optimization. Currently the local derivation free algorithms form NLOPT pack were tested. The methods: `:LN_NELDERMEAD, :LN_COBYLA, :LN_PRAXIS` show good results. Methods: `:LN_BOBYQA, :LN_SBPLX, :LN_NEWUOA` is not recommended.
- `kwargs...`: *experimental*. other keyword arguments passed to optimization method.
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
    loss_tol::Float64 = 1e-3, # i do not know how to use it
    local_alg::Symbol = :LN_NELDERMEAD,
    kwargs... # options for local fitter
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
