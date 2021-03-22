"""
    logit10(x::Real)
Function transforming interval [0,1] to [-Inf, Inf] using logit transformation.
"""
function logit10(x::Real)
    log10(x / (1.0 - x))
end

"""
    logistic10(x::Real)
Function transforming interval [-Inf, Inf] to [0,1] using logistic transformation.
Inverse function for [`logit10`](@ref).
"""
function logistic10(x::Real)
    exp10(x) / (exp10(x) + 1.0)
end

"""
    scaling(x::Real, scale::Symbol = :direct)
Transforms values from specific scale to range [-Inf, Inf] based on option.

## Return
Transformed value.

## Arguments
* `x`: input value.
* `scale`: transformation type: `:direct, :log, :logit`.
"""
function scaling(x::Real, scale::Symbol = :direct)
    if scale == :direct
        return x
    elseif scale == :log
        return log10(x)
    elseif scale == :logit
        return logit10(x)
    else
        throw(DomainError(scale, "scale type is not supported"))
    end
end

"""
    unscaling(x::Real, scale::Symbol = :direct)
Transforms values from [-Inf, Inf] to specific scale based on option. Inverse function
for [`scaling`](@ref).

## Return
Transformed value.

## Arguments
* `x`: input value.
* `scale`: transformation type: `:direct, :log, :logit`.
"""
function unscaling(x::Real, scale::Symbol = :direct)
    if scale == :direct
        return x
    elseif scale == :log
        return exp10(x)
    elseif scale == :logit
        return logistic10(x)
    else
        throw(DomainError(scale, "scale type is not supported"))
    end
end

unscaling(x::T, scale::Symbol) where T <: ForwardDiff.Dual = x

scaling(x_tup::Tuple{Float64,Float64}, scale::Symbol = :direct) = scaling.(x_tup, scale)

unscaling(x_tup::Tuple{Float64,Float64}, scale::Symbol = :direct) = unscaling.(x_tup, scale)

scaling(::Nothing, ::Symbol) = nothing

unscaling(::Nothing, ::Symbol) = nothing

"""
    function get_endpoint(
        theta_init::Vector{Float64},
        theta_num::Int,
        loss_func::Function,
        method::Symbol,
        direction::Symbol = :right;

        loss_crit::Float64 = 0.0,
        scale::Vector{Symbol} = fill(:direct, length(theta_init)),
        theta_bounds::Vector{Tuple{Float64,Float64}} = unscaling.(
            fill((-Inf, Inf), length(theta_init)),
            scale
            ),
        scan_bound::Float64 = unscaling(
            (direction==:left) ? -9.0 : 9.0,
            scale[theta_num]
            ),
        scan_tol::Float64 = 1e-3,
        loss_tol::Float64 = 1e-3,
        local_alg::Symbol = :LN_NELDERMEAD,
        kwargs...
        )

Calculates right or left endpoint of CI for parameter component. It is a wripper
of `get_right_endpoint` functions for selection of direction and using different
transformations for faster optimization.

## Return
[`EndPoint`](@ref) object storing confidence endpoint and profile points found on fly.

## Arguments
- `theta_init`: starting values of parameter vector ``\\theta``. The starting values is not necessary to be the optimum values for `loss_func` but it the value of `loss_func` must be lower than `loss_crit`.
- `theta_num`: number ``n`` of vector component to compute confidence interval ``\\theta^n``.
- `loss_func`: loss function ``\\Lambda\\left(\\theta\\right)`` the profile of which is analyzed. Usually we use log-likelihood for profile analysis in form ``\\Lambda( \\theta ) = - 2 ln\\left( L(\\theta) \\right)``.
- `method`: computational method to evaluate interval endpoint. Currently the following methods are implemented: `:CICO_ONE_PASS`, `:LIN_EXTRAPOL`, `:QUADR_EXTRAPOL`.
- `direction`: `:right` or `:left` endpoint to estimate.

## Keyword arguments
see [`get_interval`](@ref)

"""
function get_endpoint(
    theta_init::Vector{Float64},
    theta_num::Int,
    loss_func::Function,
    method::Symbol,
    direction::Symbol = :right;

    loss_crit::Float64 = 0.0,
    # :direct, :log, :logit
    scale::Vector{Symbol} = fill(:direct, length(theta_init)),
    theta_bounds::Vector{Tuple{Float64,Float64}} = unscaling.(
        fill((-Inf, Inf), length(theta_init)),
        scale
        ),
    scan_bound::Float64 = unscaling(
        (direction==:left) ? -9.0 : 9.0,
        scale[theta_num]
        ),
    scan_tol::Float64 = 1e-3,
    loss_tol::Float64 = 1e-3,
    local_alg::Symbol = :LN_NELDERMEAD,
    kwargs... # other options for get_right_endpoint
    )
    isLeft = direction == :left

    # checking arguments
    # theta_bound[1] < theta_init < theta_bound[2]
    theta_init_outside_theta_bounds = .! [theta_bounds[i][1] < theta_init[i] < theta_bounds[i][2] for i in 1:length(theta_init)]
    if any(theta_init_outside_theta_bounds)
        throw(ArgumentError("theta_init is outside theta_bound: $(findall(theta_init_outside_theta_bounds))"))
    end
    # scan_bound should be within theta_bounds
    !(theta_bounds[theta_num][1] < scan_bound < theta_bounds[theta_num][2]) &&
        throw(ArgumentError("scan_bound are outside of the theta_bounds $(theta_bounds[theta_num])"))
    # theta_init should be within scan_bound
    if (theta_init[theta_num] >= scan_bound && !isLeft) || (theta_init[theta_num] <= scan_bound && isLeft)
        throw(ArgumentError("init values are outside of the scan_bound $scan_bound"))
    end
    # 0 <= theta_bound[1] for :log
    less_than_zero_theta_bounds = (scale .== :log) .& [theta_bounds[i][1] < 0 for i in 1:length(theta_init)]
    if any(less_than_zero_theta_bounds)
        throw(ArgumentError(":log scaled theta_bound min is negative: $(findall(less_than_zero_theta_bounds))"))
    end
    # 0 <= theta_bounds <= 1 for :logit
    less_than_zero_theta_bounds = (scale .== :logit) .& [theta_bounds[i][1] < 0 || theta_bounds[i][2] > 1 for i in 1:length(theta_init)]
    if any(less_than_zero_theta_bounds)
        throw(ArgumentError(":logit scaled theta_bound min is outside range [0,1]: $(findall(less_than_zero_theta_bounds))"))
    end
    # loss_func(theta_init) < loss_crit
    !(loss_func(theta_init) < loss_crit) &&
        throw(ArgumentError("Check theta_init and loss_crit: loss_func(theta_init) should be < loss_crit"))

    # set counter in the scope
    counter::Int = 0
    # set supreme, maximal or minimal value of scanned parameter inside critical
    supreme_gd = nothing

    # transforming
    theta_init_gd = scaling.(theta_init, scale)
    if isLeft theta_init_gd[theta_num] *= -1 end # change direction
    function loss_func_gd(theta_gd)
        theta_g = copy(theta_gd)
        if isLeft theta_g[theta_num] *= -1 end # change direction
        theta = unscaling.(theta_g, scale)
        # calculate function
        loss_norm = loss_func(theta) - loss_crit

        # update counter
        counter += 1
        # update supreme ?
        update_supreme = (loss_norm < 0.) &&
            (typeof(supreme_gd)==Nothing || (theta_gd[theta_num] > supreme_gd))
        if update_supreme
            supreme_gd = theta_gd[theta_num]
        end

        return loss_norm
    end
    theta_bounds_gd = scaling.(theta_bounds, scale)
    if isLeft theta_bounds_gd[theta_num] = (-1*theta_bounds_gd[theta_num][2], -1*theta_bounds_gd[theta_num][1]) end # change direction
    scan_bound_gd = scaling(scan_bound, scale[theta_num])
    if isLeft scan_bound_gd *= -1 end # change direction
    
    # calculate endpoint using base method
    (optf_gd, pp_gd, status) = get_right_endpoint(
        theta_init_gd,
        theta_num,
        loss_func_gd,
        Val(method);
        theta_bounds = theta_bounds_gd,
        scan_bound = scan_bound_gd,
        scan_tol = scan_tol,
        loss_tol = loss_tol,
        local_alg = local_alg,
        kwargs... # options for local fitter
    )

    # transforming back
    if (isLeft && typeof(optf_gd)!==Nothing) optf_gd *= -1 end # change direction
    optf = unscaling(optf_gd, scale[theta_num])
    temp_fun = (pp::ProfilePoint) -> begin
        if isLeft pp.params[theta_num] *= -1 end # change direction
        ProfilePoint(
            unscaling(pp.params[theta_num], scale[theta_num]),
            pp.loss + loss_crit,
            unscaling.(pp.params, scale),
            pp.ret,
            pp.counter
        )
    end
    pps = [ temp_fun(pp_gd[i]) for i in 1:length(pp_gd) ]
    # transforming supreme back
    if (isLeft && typeof(supreme_gd)!==Nothing) supreme_gd *= -1 end # change direction
    supreme = unscaling(supreme_gd, scale[theta_num])

    EndPoint(optf, pps, status, direction, counter, supreme)
end


"""
    function get_endpoint(
        theta_init::Vector{Float64},
        scan_func::Function,
        loss_func::Function,
        method::Symbol,
        direction::Symbol = :right;

        loss_crit::Float64 = 0.0,
        scale::Vector{Symbol} = fill(:direct, length(theta_init)),
        theta_bounds::Vector{Tuple{Float64,Float64}} = unscaling.(
            fill((-Inf, Inf), length(theta_init)),
            scale
            ),
        scan_bound::Float64 = unscaling(
            (direction==:left) ? -9.0 : 9.0,
            :direct
            ),
        scan_tol::Float64 = 1e-3,
        loss_tol::Float64 = 1e-3,
        local_alg::Symbol = :LN_NELDERMEAD,
        kwargs...
        )

Calculates right or left endpoint of CI for parameter component. It is a wripper
of `get_right_endpoint` functions for selection of direction and using different
transformations for faster optimization.

## Return
[`EndPoint`](@ref) object storing confidence endpoint and profile points found on fly.

## Arguments
- `theta_init`: starting values of parameter vector ``\\theta``. The starting values is not necessary to be the optimum values for `loss_func` but it the value of `loss_func` must be lower than `loss_crit`.
- `scan_func`: scan function of parameters vector.
- `loss_func`: loss function ``\\Lambda\\left(\\theta\\right)`` the profile of which is analyzed. Usually we use log-likelihood for profile analysis in form ``\\Lambda( \\theta ) = - 2 ln\\left( L(\\theta) \\right)``.
- `method`: computational method to evaluate interval endpoint. Currently the following methods are implemented: `:CICO_ONE_PASS`, `:LIN_EXTRAPOL`, `:QUADR_EXTRAPOL`.
- `direction`: `:right` or `:left` endpoint to estimate.

## Keyword arguments
see [`get_interval`](@ref)

"""
function get_endpoint(
    theta_init::Vector{Float64},
    scan_func::Function,
    loss_func::Function,
    method::Symbol,
    direction::Symbol = :right;

    loss_crit::Float64 = 0.0,
    # :direct, :log, :logit
    scale::Vector{Symbol} = fill(:direct, length(theta_init)),
    theta_bounds::Vector{Tuple{Float64,Float64}} = unscaling.(
        fill((-Inf, Inf), length(theta_init)),
        scale
        ),
    scan_bound::Float64 = (direction==:left) ? -1e9 : 1e9, # log scan bound is not implemented
    scan_tol::Float64 = 1e-3,
    loss_tol::Float64 = 1e-3,
    local_alg::Symbol = :LN_NELDERMEAD,
    kwargs... # other options for get_right_endpoint
    )
    isLeft = direction == :left

    # checking arguments
    # theta_bound[1] < theta_init < theta_bound[2]
    theta_init_outside_theta_bounds = .! [theta_bounds[i][1] < theta_init[i] < theta_bounds[i][2] for i in 1:length(theta_init)]
    if any(theta_init_outside_theta_bounds)
        throw(ArgumentError("theta_init is outside theta_bound: $(findall(theta_init_outside_theta_bounds))"))
    end

    # 0 <= theta_bound[1] for :log
    less_than_zero_theta_bounds = (scale .== :log) .& [theta_bounds[i][1] < 0 for i in 1:length(theta_init)]
    if any(less_than_zero_theta_bounds)
        throw(ArgumentError(":log scaled theta_bound min is negative: $(findall(less_than_zero_theta_bounds))"))
    end
    # 0 <= theta_bounds <= 1 for :logit
    less_than_zero_theta_bounds = (scale .== :logit) .& [theta_bounds[i][1] < 0 || theta_bounds[i][2] > 1 for i in 1:length(theta_init)]
    if any(less_than_zero_theta_bounds)
        throw(ArgumentError(":logit scaled theta_bound min is outside range [0,1]: $(findall(less_than_zero_theta_bounds))"))
    end
    # loss_func(theta_init) < loss_crit
    !(loss_func(theta_init) < loss_crit) &&
        throw(ArgumentError("Check theta_init and loss_crit: loss_func(theta_init) should be < loss_crit"))

    # set counter in the scope
    counter::Int = 0
    # set supreme, maximal or minimal value of scanned parameter inside critical
    supreme_gd = nothing

    # transforming
    theta_init_gd = scaling.(theta_init, scale)
    
    function scan_func_gd(theta_gd)
        theta_g = copy(theta_gd)
        theta = unscaling.(theta_g, scale)
        scan_val = scan_func(theta)
        scan_val_gd = isLeft ? (-1)*scan_val : scan_val

        return scan_val_gd
    end

    function loss_func_gd(theta_gd)
        theta_g = copy(theta_gd)
        theta = unscaling.(theta_g, scale)
        # calculate function
        loss_norm = loss_func(theta) - loss_crit

        # update counter
        counter += 1
        # update supreme ?
        #scan_val = scan_func(theta)
        #scan_val_gd = isLeft ? (-1)*scan_val : scan_val
        #update_supreme = (loss_norm < 0.) &&
        #    (typeof(supreme_gd)==Nothing || (scan_val_gd > supreme_gd))
        #if update_supreme
        #    supreme_gd = scan_val_gd
        #end

        return loss_norm
    end
    theta_bounds_gd = scaling.(theta_bounds, scale)

    # TODO: transformed by scan_scale: scaling(scan_bound, scan_scale)
    scan_bound_gd = scan_bound
    if isLeft scan_bound_gd *= -1 end # change direction

    # calculate endpoint using base method
    (optf_gd, pp_gd, status) = get_right_endpoint(
        theta_init_gd,
        scan_func_gd,
        loss_func_gd,
        Val(method);
        theta_bounds = theta_bounds_gd,
        scan_bound = scan_bound_gd,
        scan_tol = scan_tol,
        loss_tol = loss_tol,
        local_alg = local_alg,
        kwargs... # options for local fitter
    )

    # transforming back
    temp_fun = (pp::ProfilePoint) -> begin
        value = isLeft ? (-1)*pp.value : pp.value # TODO: scan_scale
        ProfilePoint(
            value,
            pp.loss + loss_crit,
            unscaling.(pp.params, scale),
            pp.ret,
            pp.counter
        )
    end
    pps = [ temp_fun(pp_gd[i]) for i in 1:length(pp_gd) ]
    # transforming supreme back
    if (isLeft && typeof(supreme_gd)!==Nothing) supreme_gd *= -1 end # change direction
    supreme = unscaling(supreme_gd, :direct) # scan_scale ??

    if (isLeft && typeof(optf_gd)!==Nothing) optf_gd *= -1 end # change direction
    # optf = unscaling(optf_gd, scan_scale)
    optf = optf_gd

    EndPoint(optf, pps, status, direction, counter, supreme)
end

"""
    function get_right_endpoint(
        theta_init::Vector{Float64},
        theta_num::Int,
        loss_func::Function,
        method::Val{:CICO_ONE_PASS};

        theta_bounds::Vector{Tuple{Float64,Float64}} = fill(
            (-Inf, Inf), length(theta_init)
            ),
        scan_bound::Float64 = 9.0,
        scan_tol::Float64 = 1e-3,
        loss_tol::Float64 = 0.,
        local_alg::Symbol = :LN_NELDERMEAD,
        kwargs...
        )
Interface for current and future methods for endpoint estimation.

## Return
Tuple of three values:
- Right end point value: `::Float64`.
- Profile points estimated on fly: `::Array{ ProfilePoint, 1}`, see [`ProfilePoint`](@ref).
- Status of sulution: `::Symbol`. One of values: `:BORDER_FOUND_BY_SCAN_TOL`, `:SCAN_BOUND_REACHED`.

## Arguments
- `theta_init`: starting values of parameter vector ``\\theta``. The starting values is not necessary to be the optimum values for `loss_func` but it the value of `loss_func` must be lower than `loss_crit`.
- `theta_num`: number ``n`` of vector component to compute confidence interval ``\\theta^n``.
- `loss_func`: loss function the profile of which is analyzed, see [`get_interval`](@ref). In this `function loss` crit is always equal 0 for code simplification.
- `method`: this value is always fixed. Implemented methods are: `Val{:CICO_ONE_PASS}`. It is implemented for easy switching between different implemented and future methods.

## Keyword arguments
- `theta_bound`: vector of bounds for each component in format `(left_bound, right_bound)`. This bounds define the ranges for possible parameter values.
- `scan_bound`: right scan bound for `theta_num` component. It must be within the `theta_bounds` for the scanned component.
- `scan_tol`: Absolute tolerance of scanned component (stop criterion).
- `loss_tol`: Absolute tolerance of `loss_func` at `loss_crit` (stop criterion). *Restriction*. Currently is not effective for `:CICO_ONE_PASS` methods because of limitation in `LN_AUGLAG` interface.
- `local_alg`: algorithm of optimization. Currently the local derivation free algorithms form NLOPT pack were tested. The methods: `:LN_NELDERMEAD, :LN_COBYLA, :LN_PRAXIS` show good results. Methods: `:LN_BOBYQA, :LN_SBPLX, :LN_NEWUOA` is not recommended.
- `kwargs...`: the additional keyword arguments passed to `get_right_endpoint` for specific `method`.
"""
function get_right_endpoint
end
