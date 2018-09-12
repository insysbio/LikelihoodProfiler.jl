# Pkg.add("NLopt")

using NLopt, Parameters

garmonize(x::Float64, logscale::Bool) = logscale ? log10(x) : x
garmonize(b::Vector{Float64}, logscale::Bool) = logscale ? log10.(b) : b

ungarmonize(x::Float64, logscale::Bool) = logscale ? exp10(x) : x
ungarmonize(x::Vector{Float64}, logscale::Bool) = logscale ? exp10.(x) : x

"Structure storing one point from profile"
struct ProfilePoint
    loss::Float64
    params::Array{Float64, 1}
end

"Structure storing input data for parameter interval calculation"
struct ParamInput
    init_params::Vector{Float64} # initial parameters vector
    id::Int64 # id of the parameter for analysis
    loss_crit::Float64 # loss function maximum value, "identifiability level"
    loss_func::Function # loss function

    logscale::Vector{Bool} # bool vector length(init_params) where true - log scale / false - direct scale
    scan_bound::Vector{Float64} # search bounds for id parameter (default - [1e-9,1e9])
    bounds::Vector{Vector{Float64}} # bound constraints for all parameters except id
    local_alg::Symbol # local fitting algorithm (default - :LN_NELDERMEAD)
    max_iter::Int64 # maximum function evaluations
    ptol::Float64 # fitting tolerance for local optimizer (default - 1e-3)
    losstol::Float64 # constraints tolerance
end

"Structure storing result of parameter interval calculation"
struct ParamInterval
    intervals::Array{Float64, 1} # result of interval calculation. If open intervals than undefined
    ret_codes::Array{Symbol, 1} # returned result: :BOUNDS_REACHED if cannot calculate, :FTOL_REACHED if everything ok
    count_evals::Array{Int64, 1} # count of loss_func calls
    loss_final::Array{Float64, 1} # value of loss_func calculated on intervals or scan_bound

	method::Symbol # method of interval calculation: :ONE_PASS (our method)
	loss_crit::Float64 # critical level of loss function
	scan_bound::Array{Float64,1} # parameteer bound for scan
	local_alg::Symbol # algorythm of fitting, now tested on :LN_NELDERMEAD
	ptol::Float64 # required tolerance for interval estimation
    losstol::Float64 # required tolerance for loss_function at intervals

    profile_buffer::Array{ProfilePoint, 1} # storage for true points of loss_func profile
end

"""
# Input:
    init_params - initial parameters vector
    id - id of the parameter for analysis
    loss_crit - loss function maximum value, "identifiability level"
    loss_func - loss function
    logscale_all - set logscale for all parameters to true / false
    logscale - bool vector length(init_params) where true - log scale / false - direct scale
    scan_bound - search bounds for id parameter (default - [1e-9,1e9])
    fit_alg - fitting algorithm (default - :LN_AUGLAG)
    local_alg - local fitting algorithm (default - :LN_NELDERMEAD)
    bounds - bound constraints for all parameters except id
    max_iter - maximum function evaluations
    ptol - fitting tolerance for local optimizer (default - 1e-3)
    losstol - constraints tolerance

# Return:
    confidence intervals evaluation:
    (interval, termination reason, numer of evaluations, loss value)
"""
function params_intervals(
    init_params::Vector{Float64},
    id::Int64,
    loss_crit::Float64,
    loss_func::Function; # (Array{Float64,1})

    method::Symbol = :ONE_PASS,
    logscale_all::Bool = false,
    logscale::Vector{Bool} = fill(logscale_all, length(init_params)),
    scan_bound::Vector{Float64} = ungarmonize.(
        [-9., 9.],
        logscale[id]
    ),
    # fit_alg::Symbol = :LN_AUGLAG,
    bounds::Vector{Vector{Float64}} = ungarmonize.(
        fill([-Inf, Inf], length(init_params)),
        logscale
    ),
    local_alg::Symbol = :LN_NELDERMEAD,
    max_iter::Int64 = 100000,
    ptol::Float64 = 1e-3,
    losstol::Float64 = 1e-3, # tol_const::Float64 = 1e-3
    # ftol_glob::Float64 = 0.
)

    # Checking arguments
    # init_params
    !(loss_func(init_params) < loss_crit) &&
        throw(ArgumentError("Check init_params and loss_crit: loss_func(init_params) should be < loss_crit"))
    # scan bounds should be within bounds
    !(bounds[id][1] < scan_bound[1] < scan_bound[2] < bounds[id][2]) &&
        throw(ArgumentError("scan bounds are outside of the bounds $bound[id]"))
    # init_params should be within scan_bound
    !(scan_bound[1] < init_params[id] < scan_bound[2]) &&
        throw(ArgumentError("init values are outside of the scan_bound $scan_bound"))


    # Input
    input = ParamInput(
        init_params,
        id,
        loss_crit,
        loss_func,

        logscale,
        scan_bound,
        bounds,
        local_alg,
        max_iter,
        ptol,
        losstol
    )

    # Output
    result = ParamInterval(
        Vector{Float64}(2),
        Vector{Symbol}(2),
        Vector{Int64}(2),
        Vector{Float64}(2),

    	method,
    	loss_crit,
    	scan_bound,
    	local_alg,
    	ptol,
        losstol,

        []
    )

    interval_calc(input, Val(method), result)
end
