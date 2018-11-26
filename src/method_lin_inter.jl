import LikelihoodProfiler: get_right_endpoint

function get_right_endpoint(
    theta_init::Vector{Float64}, # initial point of parameters
    theta_num::Int, # number of parameter to scan
    loss_func::Function, # lambda(theta) - labmbda_min - delta_lambda
    method::Val{:LIN_INTER}; # function works only for method LIN_INTER

    theta_bounds::Vector{Tuple{Float64,Float64}} = fill(
        (-Inf, Inf), length(theta_init)
        ),
    scan_bound::Float64 = 9.0,
    scan_tol::Float64 = 1e-3,
    loss_tol::Float64 = 1e-3,
    # method args
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

    pps = ProfilePoint[]
    i = 1

    # first iteration
    current_x = theta_init[theta_num]
    current_point = prof(current_x)
    push!(pps, current_point)
    # next step
    interpolated_x = Inf
    next_x = minimum([current_x+scan_hmax, interpolated_x, scan_bound])
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
        elseif isapprox(current_x, previous_x, atol = scan_tol) && sign(current_x)*sign(current_x)<=0. # TODO: bad results are possible
            return (current_x, pps, :BORDER_FOUND_BY_SCAN_TOL) # break
        end
        # next step
        if (current_point.loss - previous_point.loss) / (current_x - previous_x) <= 0
            interpolated_x = Inf
        else
            interpolated_x = previous_x - (current_x - previous_x) * previous_point.loss / (current_point.loss - previous_point.loss)
        end
        next_x = minimum([current_x+scan_hmax, interpolated_x, scan_bound])
        (previous_point, previous_x, current_x) = (current_point, current_x, next_x)
    end

    return nothing
end
