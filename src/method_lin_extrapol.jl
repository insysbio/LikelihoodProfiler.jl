import LikelihoodProfiler: get_right_endpoint

function get_right_endpoint(
    theta_init::Vector{Float64}, # initial point of parameters
    theta_num::Int, # number of parameter to scan
    loss_func::Function, # lambda(theta) - labmbda_min - delta_lambda
    method::Val{:LIN_EXTRAPOL}; # function works only for method LIN_INTER

    theta_bounds::Vector{Tuple{Float64,Float64}} = fill(
        (-Inf, Inf), length(theta_init)
        ),
    scan_bound::Float64 = 9.0,
    scan_tol::Float64 = 1e-3,
    loss_tol::Float64 = 1e-3,
    # method args
    scan_hini = 1.,
    scan_hmax = Inf,
    # local alg args
    local_alg::Symbol = :LN_NELDERMEAD,
    max_iter::Int = 10^5,
    ftol_abs::Float64 = 1e-3,
    kwargs... # options for local fitter
    )
    # dim of the theta vector
    n_theta = length(theta_init)

    prof = profile(
        theta_init,
        theta_num,
        loss_func;
        theta_bounds = theta_bounds,
        local_alg = local_alg,
        ftol_abs = loss_tol,
        maxeval = max_iter
    )

    # empty container
    pps = ProfilePoint[]

    # first iteration
    current_x = theta_init[theta_num]
    current_point = prof(current_x)
    push!(pps, current_point)

    i = 1
    # next step
    extrapol_x = current_x + scan_hini
    next_x = minimum([current_x+scan_hmax, extrapol_x])
    (previous_point, previous_x, current_x) = (current_point, current_x, next_x)

    # other iterations
    while true
        current_point = prof(current_x)
        push!(pps, current_point)

        i += 1
        if i > 100
            return (scan_bound, pps, :MAX_ITER_REACHED) # break
        elseif current_x >= scan_bound && current_point.loss < 0.
            return (scan_bound, pps, :SCAN_BOUND_REACHED) # break
        elseif isapprox(current_point.loss, 0., atol = loss_tol)
            return (current_x, pps, :BORDER_FOUND_BY_LOSS_TOL) # break
        elseif isapprox(current_x, previous_x, atol = scan_tol)
            return (current_x, pps, :BORDER_FOUND_BY_SCAN_TOL) # break
        end

        # next step
        if (current_point.loss - previous_point.loss) / (current_x - previous_x) <= 0
            extrapol_x = current_x + scan_hini
        else
            extrapol_x = previous_x - (current_x - previous_x) * previous_point.loss / (current_point.loss - previous_point.loss)
        end
        next_x = minimum([current_x+scan_hmax, extrapol_x])
        (previous_point, previous_x, current_x) = (current_point, current_x, next_x)
    end
end
