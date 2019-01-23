
using LinearAlgebra

function get_right_endpoint(
    theta_init::Vector{Float64}, # initial point of parameters
    theta_num::Int, # number of parameter to scan
    loss_func::Function, # lambda(theta) - labmbda_min - delta_lambda
    method::Val{:QUADR_EXTRAPOL}; # function works only for method LIN_INTER

    theta_bounds::Vector{Tuple{Float64,Float64}} = fill(
        (-Inf, Inf), length(theta_init)
        ),
    scan_bound::Float64 = 9.0,
    scan_tol::Float64 = 1e-3,
    loss_tol::Float64 = 0., # 1e-3,
    # method args
    scan_hini = 1.,
    scan_hmax = Inf,
    # local alg args
    local_alg::Symbol = :LN_NELDERMEAD,
    max_iter::Int = 10^5,
    ftol_abs::Float64 = 1e-3,
    kwargs... # options for local fitter
    )
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
    x_3 = theta_init[theta_num]
    theta_init_3 = theta_init
    iteration_count = 0

    # other iterations
    while true
        # preparation
        global x_2 # not initialized for the first iteration
        global x_1 # not initialized for the first iteration
        global point_2 # not initialized for the first iteration
        global point_1 # not initialized for the first iteration
        iteration_count += 1 # to understand if this is a first iteration

        # get profile point
        point_3 = prof(
            x_3;
            theta_init_i = theta_init_3, # hypothetically this makes optimization more effective
            maxeval = max_iter - accum_counter # how many calls left
            )
        push!(pps, point_3)
        accum_counter += point_3.counter # update counter

        if x_3 >= scan_bound && point_3.loss < 0.
            return (nothing, pps, :SCAN_BOUND_REACHED) # break
        elseif isapprox(point_3.loss, 0., atol = loss_tol)
            return (x_3, pps, :BORDER_FOUND_BY_LOSS_TOL) # break
        # no checking for the first and second iterations
        # elseif iteration_count>2 && isapprox(x_3, x_2, atol = scan_tol)
        elseif iteration_count>1 && isapprox((x_3 - x_2) * point_3.loss / (point_3.loss - point_2.loss), 0., atol = scan_tol)
            return (x_3, pps, :BORDER_FOUND_BY_SCAN_TOL) # break
        elseif point_3.ret == :MAXEVAL_REACHED
            return (nothing, pps, :MAX_ITER_REACHED) # break
        end

        # next step
        if iteration_count==1
            x_4_extrapol = x_3 + scan_hini
            x_4 = minimum([x_3+scan_hmax, x_4_extrapol, theta_bounds[theta_num][2]])
            (point_2, x_2, x_3, theta_init_3) = (point_3, x_3, x_4, point_3.params)
        elseif iteration_count==2
            D = [
                x_3^2 x_3 1;
                x_2^2 x_2 1;
                2*x_2 1 0
                ]
            Da = [
                point_3.loss x_3 1;
                point_2.loss x_2 1;
                0 1 0
                ]
            Db = [
                x_3^2 point_3.loss 1;
                x_2^2 point_2.loss 1;
                2*x_2 0 0
                ]
            Dc = [
                x_3^2 x_3 point_3.loss;
                x_2^2 x_2 point_2.loss;
                2*x_2 1 0
                ]
            a = det(Da) / det(D)
            b = det(Db) / det(D)
            c = det(Dc) / det(D)
            if a <= 0.
                x_4_extrapol = x_3 + scan_hini
            else
                x_4_extrapol = (-b + sqrt(b^2-4*a*c))/2/a
            end
            x_4 = minimum([x_3+scan_hmax, x_4_extrapol, theta_bounds[theta_num][2]])
            (x_1, x_2, x_3, point_1, point_2, theta_init_3) = (x_2, x_3, x_4, point_2, point_3, point_3.params)
        else
            D = [
                x_3^2 x_3 1;
                x_2^2 x_2 1;
                x_1^2 x_1 1
                ]
            Da = [
                point_3.loss x_3 1;
                point_2.loss x_2 1;
                point_1.loss x_1 1
                ]
            Db = [
                x_3^2 point_3.loss 1;
                x_2^2 point_2.loss 1;
                x_1^2 point_1.loss 1
                ]
            Dc = [
                x_3^2 x_3 point_3.loss;
                x_2^2 x_2 point_2.loss;
                x_1^2 x_1 point_1.loss
                ]
            a = det(Da) / det(D)
            b = det(Db) / det(D)
            c = det(Dc) / det(D)
            if a <= 0. || b^2-4*a*c < 0.
                x_4_extrapol = x_3 + scan_hini
            else
                x_4_extrapol = (-b + sqrt(b^2-4*a*c))/2/a
            end
            x_4 = minimum([x_3+scan_hmax, x_4_extrapol, theta_bounds[theta_num][2]])
            (x_1, x_2, x_3, point_1, point_2, theta_init_3) = (x_2, x_3, x_4, point_2, point_3, point_3.params)
        end
    end
end
