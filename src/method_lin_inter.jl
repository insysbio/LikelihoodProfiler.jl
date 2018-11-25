
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

    return nothing
end
