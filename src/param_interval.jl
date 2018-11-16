struct ParamIntervalInput
    theta_init::Vector{Float64} # initial parameters vector
    theta_num::Int64 # id of the parameter for analysis
    loss_func::Function # loss function
    loss_crit::Float64 # loss function maximum value, "identifiability level"
    scale::Vector{Symbol}
    theta_bounds::Vector{Tuple{Float64, Float64}} # search bounds for id parameter
    scan_bounds::Tuple{Float64,Float64}
    scan_tol::Float64 # fitting tolerance for local optimizer (default - 1e-3)
    loss_tol::Float64 # constraints tolerance
    local_alg::Symbol # local fitting algorithm (default - :LN_NELDERMEAD)
end

struct ParamInterval
    input::ParamIntervalInput
    method::Symbol # method :CICO_ONE_PASS, :D2D_PLE
    result::Tuple{EndPoint, EndPoint}
end

function param_interval(
    theta_init::Vector{Float64},
    theta_num::Int64,
    loss_func::Function,
    method::Symbol;

    loss_crit::Float64 = 0.0,
    scale::Vector{Symbol} = fill(:direct, length(theta_init)),
    theta_bounds::Vector{Tuple{Float64,Float64}} = ungarm.(
        fill((-Inf, Inf), length(theta_init)),
        scale
        ),
    scan_bounds::Tuple{Float64,Float64} = ungarm.(
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
        local_alg
        )

    ParamInterval(
        input,
        method,
        Tuple(endpoints)
        )
end
