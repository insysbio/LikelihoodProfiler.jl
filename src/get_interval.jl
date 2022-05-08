abstract type AbstractIntervalInput end

"""
    struct ParamIntervalInput <: AbstractIntervalInput
        theta_init::Vector{Float64} # initial parameters vector
        theta_num::Int              # number of the parameter for identifiability analysis
        loss_func::Function         # loss function
        method::Symbol
        options::Any
    end

Structure storing input data for parameter identification
"""
struct ParamIntervalInput <: AbstractIntervalInput
    theta_init::Vector{Float64} # initial parameters vector
    theta_num::Int # number of the parameter for analysis
    loss_func::Function # loss function
    method::Symbol
    options::Dict{Symbol, Any}
end

struct PredictionIntervalInput <: AbstractIntervalInput
    theta_init::Vector{Float64}
    scan_func::Function
    loss_func::Function
    method::Symbol
    options::Dict{Symbol, Any}
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
    input::AbstractIntervalInput
    loss_init::Float64
    method::Symbol
    result::Tuple{EndPoint, EndPoint}
end

"""
    function get_interval(
        theta_init::Vector{Float64},
        theta_num::Int,
        loss_func::Function,
        method::Symbol;

        scale::Vector{Symbol} = fill(:direct, length(theta_init)),
        scan_bounds::Tuple{Float64,Float64} = unscaling.(
            (-9.0, 9.0),
            scale[theta_num]
            ),
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
- `scale`: vector of scale transformations for each parameters' component. Possible values: `:direct` (`:lin`), `:log`, `:logit`. This option can speed up the optimization, especially for wide `theta_bounds`. The default value is `:direct` (no transformation) for all parameters.
- `scan_bounds`: scan bounds tuple for `theta_num` parameter. Should be within the `theta_bounds` for `theta_num` parameter. Default is `(-9.,9.)` for `:direct` scales and `(1e-9, 1e+9)` for `:log`.
- `kwargs...`: the additional arguments passed to [`get_endpoint`](@ref)
"""
function get_interval(
    theta_init::Vector{Float64},
    theta_num::Int,
    loss_func::Function,
    method::Symbol;

    scale::Vector{Symbol} = fill(:direct, length(theta_init)),
    scan_bounds::Tuple{Float64,Float64} = unscaling.(
        (-9.0, 9.0),
        scale[theta_num]
        ),
    kwargs... # other options for get_right_endpoint
)
    # both endpoints
    endpoints = [get_endpoint(
        theta_init,
        theta_num,
        loss_func,
        method,
        [:left,:right][i]; # method

        scale,
        scan_bound = scan_bounds[i],
        kwargs... # options for local fitter
        ) for i in 1:2]

    input = ParamIntervalInput(
        theta_init,
        theta_num,
        loss_func,
        method,
        Dict(:scale=>scale, :scan_bounds=>scan_bounds, kwargs...) # ? NamedTuple
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
        loss_tol::Float64 = 0.,
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
- `scale`: vector of scale transformations for each parameters' component. Possible values: `:direct` (`:lin`), `:log`, `:logit`. This option can speed up the optimization, especially for wide `theta_bounds`. The default value is `:direct` (no transformation) for all parameters.
- `scan_bounds`: scan bounds tuple for `scan_func` values. Default is `(1e-9, 1e9)` .
- `kwargs...`: the additional arguments passed to [`get_endpoint`](@ref)
"""
function get_interval(
    theta_init::Vector{Float64},
    scan_func::Function,
    loss_func::Function,
    method::Symbol;

    scale::Vector{Symbol} = fill(:direct, length(theta_init)),
    scan_bounds::Tuple{Float64,Float64} = unscaling.(
        (-9.0, 9.0),
        scale[theta_num]
        ),
    kwargs... # other options for get_right_endpoint
)
    # both endpoints
    endpoints = [get_endpoint(
        theta_init,
        scan_func,
        loss_func,
        method,
        [:left,:right][i]; # method

        scale,
        scan_bound = scan_bounds[i],
        kwargs...
        ) for i in 1:2]

    input = PredictionIntervalInput(
        theta_init,
        scan_func,
        loss_func,
        method,
        Dict(:scale=>scale, :scan_bounds=>scan_bounds, kwargs...) # ? NamedTuple
        )

    ParamInterval(
        input,
        loss_func(theta_init),
        method,
        Tuple(endpoints)
        )
end
