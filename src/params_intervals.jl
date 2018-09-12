# Pkg.add("NLopt")

using NLopt

garmonize(x::Float64, logscale::Bool) = logscale ? log10(x) : x
garmonize(b::Vector{Float64}, logscale::Bool) = logscale ? log10.(b) : b

ungarmonize(x::Float64, logscale::Bool) = logscale ? exp10(x) : x
ungarmonize(x::Vector{Float64}, logscale::Bool) = logscale ? exp10.(x) : x

"Structure storing one point from profile"
struct ProfilePoint
    loss::Float64
    params::Array{Float64, 1}
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
	alg_loc::Symbol # algorythm of fitting, now tested on :LN_NELDERMEAD
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
    logscale_all - ???
    logscale - bool vector length(init_params) where true - log scale / false - direct scale
    scan_bound - search bounds for id parameter (default - [1e-9,1e9])
    fit_alg - fitting algorithm (default - :LN_AUGLAG)
    local_alg - local fitting algorithm (default - :LN_NELDERMEAD)
    bounds - bound constraints for all parameters except id
    max_iter - ???
    ptol - fitting tolerance for local optimizer (default - 1e-3)

# Return:
    confidence intervals evaluation:
    (interval, termination reason, numer of evaluations, loss value)
"""
function params_intervals(
    init_params::Vector{Float64},
    id::Int64,
    loss_crit::Float64,
    loss_func::Function; # (Array{Float64,1})

    logscale_all::Bool = false,
    logscale::Vector{Bool} = fill(logscale_all, length(init_params)),
    scan_bound::Vector{Float64} = ungarmonize.(
        [-9., 9.],
        logscale[id]
    ),
    # fit_alg::Symbol = :LN_AUGLAG,
    local_alg::Symbol = :LN_NELDERMEAD,
    bounds::Vector{Vector{Float64}} = ungarmonize.(
        fill([-Inf, Inf], length(init_params)),
        logscale
    ),
    max_iter::Int64 = 100000,
    ptol::Float64 = 1e-3,
    losstol::Float64 = 1e-3, # tol_const::Float64 = 1e-3
    # ftol_glob::Float64 = 0.
)
    # set counter scope
    counter::Int64 = 0

    # Output
    result = ParamInterval(
        Vector{Float64}(2),
        Vector{Symbol}(2),
        Vector{Int64}(2),
        Vector{Float64}(2),

		:ONE_PASS,
		loss_crit,
		scan_bound,
		local_alg,
		ptol,
        losstol,

        []
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

    # garmonize init_params
    params = garmonize.(init_params, logscale)

    # garmonize bounds
    bounds_garm = garmonize.(bounds, logscale)

    # Objective function
    # optim_func(x, g) = ungarmonize(x[id], logscale[id])
    optim_func(x, g) = x[id]

    # Constraints function
    function constraints_func(x, g)
        x_initial_scale = ungarmonize.(x, logscale) # potentiation
        loss = loss_func(x_initial_scale) - loss_crit; counter += 1

        if !(scan_bound[1] <= x_initial_scale[id] <= scan_bound[2]) && loss < 0.
            throw(ForcedStop())
        else
            return loss
        end
    end

    # Confidence interval search
    for min_max in (:min, :max)
        int_id = min_max == :min ? 1 : 2
        counter = 0 # set zero counter

        (optf, optx, ret) = params_intervals_one_side(
            params,
            optim_func,
            constraints_func,
            min_max;
            bounds = bounds_garm,
            # fit_alg = fit_alg,
            local_alg = local_alg,
            ftol_loc = ptol,
            tol_const = losstol,
            max_iter = max_iter
        )

        # if bounds reached
        if ret == :FORCED_STOP
            # result.intervals[int_id] = ungarmonize(optf, logscale[id])
            result.ret_codes[int_id] = :BOUNDS_REACHED
        else
            result.intervals[int_id] = ungarmonize(optf, logscale[id])
            local params_final = ungarmonize.(optx, logscale)
            local loss_final = loss_func(params_final); counter += 1
            push!(
                result.profile_buffer,
                ProfilePoint(loss_final, params_final)
            )
            result.ret_codes[int_id] = ret
        end

        result.loss_final[int_id] = loss_func(ungarmonize.(optx, logscale)); counter += 1
        result.count_evals[int_id] = counter
    end

    result # return
end # function

"""
# Input:
    init_params - initial parameters vector
    optim_func - ???
    constraints_func - ???
    min_max - ???

    fit_alg - fitting algorithm (default - :LN_AUGLAG)
    local_alg - local fitting algorithm (default - :LN_NELDERMEAD)
    bounds - bound constraints for all parameters except id
    max_iter - ???
    ftol_loc - fitting tolerance for local optimizer (default - 1e3)

# Return:
    parameter profile plot
"""
function params_intervals_one_side(
    init_params::Vector{Float64},
    optim_func::Function,
    constraints_func::Function,
    min_max::Symbol;

    fit_alg::Symbol = :LN_AUGLAG,
    local_alg::Symbol = :LN_NELDERMEAD,
    bounds::Vector{Vector{Float64}} = fill([-Inf, Inf], length(init_params)),
    ftol_loc::Float64 = 1e-3,
    tol_const::Float64 = 1e-3, # tolerance of constraints
    # ftol_glob::Float64 = 1e-3, # tolerance of global method
    max_iter::Int64 = 100000
)
    # dim of the problem
    n_params = length(init_params)

    # optimization obj
    opt = Opt(fit_alg, n_params)

    # min or max
    if min_max == :min
        min_objective!(opt, optim_func)
    else
        max_objective!(opt, optim_func)
    end

    # bound constraints
    lb = minimum.(bounds)
    ub = maximum.(bounds)

    # if alg = Augmented Lagrangian
    if fit_alg in (:LN_AUGLAG, :LN_AUGLAG_EQ)
        # local optimizer
        local_opt = Opt(local_alg, n_params)
        # tolerances
        # ftol_abs!(opt, ftol_glob)
        ftol_abs!(local_opt, ftol_loc)
        # bound constraints
        lower_bounds!(opt, lb)
        upper_bounds!(opt, ub)
        # adding local optimizer
        local_optimizer!(opt, local_opt)
    else  # other Algorithm
        # tolerances
        # ftol_abs!(opt, ftol_glob)
        # bound constraints
        lower_bounds!(opt, lb)
        upper_bounds!(opt, ub)
    end

    # max function calls
    maxeval!(opt, max_iter)

    # inequality constraints
    # inequality_constraint!(opt, constraints_func, )
    inequality_constraint!(opt, constraints_func, tol_const)

    # return
    (optf, optx, ret) = optimize(opt, init_params)
end # function
