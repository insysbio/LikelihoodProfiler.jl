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
    # to count loss function calls inside this function, accumulation
    accum_counter::Int = 0
    # empty container
    pps = ProfilePoint[]

    prof = profile(
        theta_init,
        theta_num,
        loss_func;
        theta_bounds = theta_bounds,
        local_alg = local_alg,
        ftol_abs = loss_tol
    )

    # first iteration
    current_x = theta_init[theta_num]
    iteration_count = 0

    # other iterations
    while true
        # preparation
        global previous_x # not initialized for the first iteration
        global previous_point # not initialized for the first iteration
        iteration_count += 1

        # get profile point
        current_point = prof(
            current_x;
            maxeval = max_iter - accum_counter # how many calls left
            )
        push!(pps, current_point)
        accum_counter += current_point.counter # update counter

        if current_x >= scan_bound && current_point.loss < 0.
            return (nothing, pps, :SCAN_BOUND_REACHED) # break
        elseif isapprox(current_point.loss, 0., atol = loss_tol)
            return (current_x, pps, :BORDER_FOUND_BY_LOSS_TOL) # break
        elseif iteration_count!=1 && isapprox(current_x, previous_x, atol = scan_tol) # no checking for the first iteration
            return (current_x, pps, :BORDER_FOUND_BY_SCAN_TOL) # break
        elseif current_point.ret == :MAXEVAL_REACHED
            return (nothing, pps, :MAX_ITER_REACHED) # break
        end

        # next step
        extrapolate_next_step =
            iteration_count!=1 && # for the first iteration
            (current_point.loss - previous_point.loss) / (current_x - previous_x) > 0
        if extrapolate_next_step
            extrapol_x = previous_x - (current_x - previous_x) * previous_point.loss / (current_point.loss - previous_point.loss)
        else
            extrapol_x = current_x + scan_hini
        end
        next_x = minimum([current_x+scan_hmax, extrapol_x])
        (previous_point, previous_x, current_x) = (current_point, current_x, next_x)
    end
end
