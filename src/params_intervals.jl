# Pkg.add("NLopt")

using NLopt

garmonize(x::Float64, logscale::Bool) = logscale ? log10(x) : x
garmonize(b::Vector{Float64}, logscale::Bool) = logscale ? log10.(b) : b

ungarmonize(x::Float64, logscale::Bool) = logscale ? exp10(x) : x
ungarmonize(x::Vector{Float64}, logscale::Bool) = logscale ? exp10.(x) : x

"Structure storing result of parameter interval calculation"
struct ParamInterval
    intervals::Array{Float64}
    ret_codes::Array{Symbol}
    count_evals::Array{Int64}
    loss_final::Array{Float64}

	method::Symbol
	loss_crit::Float64
	scan_bound::Array{Float64}
	alg_loc::Symbol
	ftol_loc::Float64
    ftol_actual::Array{Float64}
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
    ftol_loc - fitting tolerance for local optimizer (default - 1e3)

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
    fit_alg::Symbol = :LN_AUGLAG,
    local_alg::Symbol = :LN_NELDERMEAD,
    bounds::Vector{Vector{Float64}} = ungarmonize.(
        fill([-Inf, Inf], length(init_params)),
        logscale
    ),
    max_iter::Int64 = 100000,
    ftol_loc::Float64 = 1e-3,
    # ftol_glob::Float64 = 0.,
    # tol_const::Float64 = 1e-3
)
    # set counter scope
    counter::Int64 = 0

    # Output
    result = ParamInterval(
        Vector{Float64}(2),
        Vector{Symbol}(2),
        Vector{Int64}(2),
        Vector{Float64}(2),

		:one,
		loss_crit,
		scan_bound,
		local_alg,
		ftol_loc,
        Vector{Float64}(2)
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
    for minmax in (:min, :max)
        int_id = minmax == :min ? 1 : 2
        counter = 0

        (optf, optx, ret) = params_intervals_one_side(
            params,
            optim_func,
            constraints_func,
            minmax;
            bounds = bounds_garm,
            fit_alg = fit_alg,
            local_alg = local_alg,
            ftol_loc = ftol_loc,
            # ftol_glob::Float64 = 1e-3, # tolerance of global method
            # tol_const::Float64 = 1e-3 # tolerance of constraints
            max_iter = max_iter
        )

        # if bounds reached
        if ret == :FORCED_STOP
            result.intervals[int_id] = minmax == :min ? scan_bound[1] : scan_bound[2]
            result.ret_codes[int_id] = :BOUNDS_REACHED
        else
            result.intervals[int_id] = ungarmonize(optf, logscale[id])
            result.ret_codes[int_id] = ret
        end

        result.loss_final[int_id] = loss_func(ungarmonize.(optx, logscale)); counter += 1

        result.ftol_actual[int_id] = abs(result.loss_final[int_id] - result.loss_crit)
        result.count_evals[int_id] = counter
    end

    result # return
end # function

"""
# Input:
    init_params - initial parameters vector
    optim_func - ???
    constraints_func - ???
    minmax - ???

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
    minmax::Symbol;

    fit_alg::Symbol = :LN_AUGLAG,
    local_alg::Symbol = :LN_NELDERMEAD,
    bounds::Vector{Vector{Float64}} = fill([-Inf, Inf], length(init_params)),
    max_iter::Int64 = 100000,
    ftol_loc::Float64 = 1e-3
    # ftol_glob::Float64 = 1e-3, # tolerance of global method
    # tol_const::Float64 = 1e-3 # tolerance of constraints
)
    # dim of the problem
    n_params = length(init_params)

    # optimization obj
    opt = Opt(fit_alg, n_params)

    # min or max
    if minmax == :min
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
    inequality_constraint!(opt, constraints_func)

    # return
    (optf, optx, ret) = optimize(opt, init_params)
end # function
