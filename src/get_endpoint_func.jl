
function get_endpoint(
    theta_init::Vector{Float64},
    #scan_func::Function, # h(theta) function for predictions or parameter
    #loss_func::Function,
    scan_loss_func::Function, # function returns tuple of scan value, loss value: (h, lambda)
    method::Symbol, # method::Val{:CICO_ONE_PASS},
    direction::Symbol = :right;

    loss_crit::Float64 = 0.0,
    scan_scale::Symbol = :direct, # scale for scan_value XXX: not works
    # :direct, :log, :logit
    scale::Vector{Symbol} = fill(:direct, length(theta_init)),
    theta_bounds::Vector{Tuple{Float64,Float64}} = unscaling.(
        fill((-Inf, Inf), length(theta_init)),
        scale
        ),
    scan_bound::Float64 = unscaling(
        (direction==:left) ? -9.0 : 9.0,
        scan_scale
        ),
    scan_tol::Float64 = 1e-3,
    loss_tol::Float64 = 1e-3, # does not work for CICO
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

    # 0 <= theta_bound for :log
    less_than_zero_theta_bounds = (scale .== :log) .& [theta_bounds[i][1] < 0 for i in 1:length(theta_init)]
    if any(less_than_zero_theta_bounds)
        throw(ArgumentError(":log scaled theta_bound min is negative: $(findall(less_than_zero_theta_bounds))"))
    end
    # 0 <= theta_bounds <= 1 for :logit
    less_than_zero_theta_bounds = (scale .== :logit) .& [theta_bounds[i][1] < 0 || theta_bounds[i][2] > 1 for i in 1:length(theta_init)]
    if any(less_than_zero_theta_bounds)
        throw(ArgumentError(":logit scaled theta_bound min is outside range [0,1]: $(findall(less_than_zero_theta_bounds))"))
    end
    # scan_initial should be within scan_bound
    # XXX: maybe it is not necessary
    #=
    scan_func_init = scan_func(theta_init)
    if (scan_func_init >= scan_bound && !isLeft) || scan_func_init <= scan_bound && isLeft)
        throw(ArgumentError("init scan_func($theta_init) is outside of the scan_bound $scan_bound"))
    end
    =#
    # loss_func(theta_init) < loss_crit
    # XXX: maybe it is not necessary
    #=
    !(loss_func(theta_init) < loss_crit) &&
        throw(ArgumentError("Check theta_init and loss_crit: loss_func(theta_init) should be < loss_crit"))
    =#
    # set counter in the scope
    counter::Int = 0
    # set supreme, maximal or minimal value of scanned parameter inside critical
    supreme_gd = nothing

    # transforming theta
    theta_init_gd = scaling.(theta_init, scale)
    # transforming scan fun
    function scan_loss_func_gd(theta_gd)
        #@show theta_gd
        theta_g = copy(theta_gd)
        theta = unscaling.(theta_g, scale)
        (scan_val, loss_val) = scan_loss_func(theta)
        scan_val_gd = isLeft ? (-1)*scan_val : scan_val

        # update counter
        counter += 1
        # update supreme
        update_supreme = (loss_val  < loss_crit) &&
            (typeof(supreme_gd) == Nothing || (scan_val_gd > supreme_gd))
        if update_supreme
            supreme_gd = scan_val_gd
        end

        (
            scan_val_gd,
            loss_val - loss_crit
        )
    end

    theta_bounds_gd = scaling.(theta_bounds, scale)

    # TODO: transformed by scan_scale: scaling(scan_bound, scan_scale)
    scan_bound_gd = scan_bound
    if isLeft scan_bound_gd *= -1 end # change direction

    # calculate endpoint using base method
    (optf_gd, pp_gd, status, iter) = get_right_endpoint(
        theta_init_gd,
        scan_loss_func_gd,
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
    supreme = unscaling(supreme_gd, scan_scale)

    if (isLeft && typeof(optf_gd)!==Nothing) optf_gd *= -1 end # change direction
    # optf = unscaling(optf_gd, scan_scale)
    optf = optf_gd

    EndPoint(optf, pps, status, direction, iter, supreme)
end
