# transforms [0,1] to [-Inf, Inf]
function logit10(x::Float64)
    log10(x / (1.0 - x))
end

# transforms [-Inf, Inf] to [0,1]
function logistic10(x::Float64)
    exp10(x) / (exp10(x) + 1.0)
end

function garm(x::Float64, scale::Symbol = :direct)
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

function ungarm(x::Float64, scale::Symbol = :direct)
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

garm(x_tup::Tuple{Float64,Float64}, scale::Symbol = :direct) = garm.(x_tup, scale)

ungarm(x_tup::Tuple{Float64,Float64}, scale::Symbol = :direct) = ungarm.(x_tup, scale)

"Structure storing one point from profile function"
struct ProfilePoint
    loss::Float64
    params::Array{Float64, 1}
    ret::Symbol
end

"End point storage"
struct EndPoint
    value::Float64
    profilePoints::Array{ProfilePoint, 1}
    status::Symbol
    direction::Symbol
    counter::Int64
end

# get left or right endpoint of CI for parameter component
function get_endpoint(
    theta_init::Vector{Float64},
    theta_num::Int64,
    loss_func::Function,
    method::Symbol,
    direction::Symbol = :right;

    loss_crit::Float64 = 0.0,
    # :direct, :log, :logit
    scale::Vector{Symbol} = fill(:direct, length(theta_init)),
    theta_bounds::Vector{Tuple{Float64,Float64}} = ungarm.(
        fill((-Inf, Inf), length(theta_init)),
        scale
        ),
    scan_bound::Float64 = ungarm(
        (direction==:left) ? -9.0 : 9.0,
        scale[theta_num]
        ),
    scan_tol::Float64 = 1e-3,
    loss_tol::Float64 = 1e-3,
    local_alg::Symbol = :LN_NELDERMEAD,
    kwargs... # options for local fitter
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
    counter::Int64 = 0

    # transforming
    theta_init_gd = garm.(theta_init, scale)
    if isLeft theta_init_gd[theta_num] *= -1 end # change direction
    function loss_func_gd(theta_gd::Vector{Float64})
        theta_g = copy(theta_gd)
        if isLeft theta_g[theta_num] *= -1 end # change direction
        theta = ungarm.(theta_g, scale)
        # update counter
        counter += 1
        # calculate function
        loss_func(theta) - loss_crit
    end
    theta_bounds_gd = garm.(theta_bounds, scale)
    if isLeft theta_bounds_gd[theta_num] = (-1*theta_bounds_gd[theta_num][2], -1*theta_bounds_gd[theta_num][1]) end # change direction
    scan_bound_gd = garm(scan_bound, scale[theta_num])
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
    if isLeft optf_gd *= -1 end # change direction
    optf = ungarm(optf_gd, scale[theta_num])
    temp_fun = (pp::ProfilePoint) -> begin
        if isLeft pp.params[theta_num] *= -1 end # change direction
        ProfilePoint(
            pp.loss + loss_crit,
            ungarm.(pp.params, scale),
            pp.ret
        )
    end
    pps = [ temp_fun(pp_gd[i]) for i in 1:length(pp_gd) ]

    EndPoint(optf, pps, status, direction, counter)
end
