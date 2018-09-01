# Pkg.add("NLopt")

using NLopt

transform_log(x::Float64, logscale::Bool) = logscale ? log10(x) : x

transform_exp(x::Float64, logscale::Bool) = logscale ? exp10(x) : x

transform_log_bounds(b::Vector{Float64}, logscale::Bool) = logscale ? Tuple(log10.(b)) : Tuple(b)

function params_intervals(
    init_params::Vector{Float64},
    id::Int64,
    maxf::Float64,
    loss_func::Function;
    logscale_all::Bool = false,
    logscale::Vector{Bool} = fill(logscale_all, length(init_params)),
    bounds_id::Vector{Float64} = [-9, 9], # [1e-9, 1e9] for log
    solver::Symbol = :NLOPT,  # currently NLOPT is the only option

    fit_alg::Symbol = :LN_AUGLAG,
    local_alg::Symbol = :LN_NELDERMEAD,
    bounds_params::Vector{Vector{Float64}} = fill(
        [-Inf, Inf], # [0, Inf] for log
        length(init_params)
    ),
    max_iter::Int64 = 100000,
    ftol_loc::Float64 = 1e-3,
    # ftol_glob::Float64 = 0.,
    # tol_const::Float64 = 1e-3
)
    # Iterations count
    global count = 0

    # Output
    intervals = Vector{Float64}(2)
    ret_codes = Vector{Symbol}(2)
    count_evals = Vector{Int64}(2)
    loss_final = Vector{Float64}(2)

    # Checks
    (loss_func(init_params) > maxf) && throw(ArgumentError("Check init_params and maxf: loss_func(init_params) should be <= maxf"))
    (init_params[id] <= minimum(bounds_id) || init_params[id] >= maximum(bounds_id)) && throw(ArgumentError("init values are outside of the bounds $bounds_id"))
    # count += 1 # because loss_func was calculated once

    # Logscale cheks - parameters
    params = transform_log.(init_params, logscale)

    # Logscale cheks - bounds
    # bounds_params[id] = bounds_id # metelkin
    bounds = transform_log_bounds.(bounds_params, logscale)

    # Objective function
    # optim_func(x, g) = transform_exp(x[id], logscale[id]) # XXX:
    optim_func(x, g) = x[id]

    # Constraints function
    function constraints_func(x, g)
        x_initial_scale = transform_exp.(x, logscale) # potentiation
        loss = loss_func(x_initial_scale) - maxf

        # global count # metelkin
        count += 1
        # print(loss, ": ")
        # println(x_initial_scale)

        if (x_initial_scale[id] <= bounds_id[1] || x_initial_scale[id] >= bounds_id[2]) && loss < 0.
            throw(ForcedStop())
        else
            return loss
        end
    end

    # Confidence interval search
    for minmax in (:min, :max)
        int_id = minmax == :min ? 1 : 2

        (optf, optx, ret) = params_intervals_one_side(
            params,
            optim_func,
            constraints_func,
            minmax;
            bounds = bounds,
            fit_alg = fit_alg,
            local_alg = local_alg,
            ftol_loc = ftol_loc,
            # ftol_glob::Float64 = 1e-3, # tolerance of global method
            # tol_const::Float64 = 1e-3 # tolerance of constraints
            max_iter = max_iter
        )

        # if bounds reached
        if ret == :FORCED_STOP
            intervals[int_id] = minmax == :min ? minimum(bounds_id) : maximum(bounds_id)
            ret_codes[int_id] = :BOUNDS_REACHED
        else
            intervals[int_id] = transform_exp(optx[id], logscale[id])
            ret_codes[int_id] = ret
        end

        loss_final[int_id] = loss_func(transform_exp.(optx, logscale))
        count += 1
        count_evals[int_id] = count

        # println("id=$id, interval[$int_id] = $(intervals[int_id]), ret_codes=$(ret_codes[int_id]), counts=$(count_evals[int_id])")
        global count = 0
    end

    # println(intervals, ", ", ret_codes, ", ", count_evals, ", ", loss_final)
    intervals, ret_codes, count_evals, loss_final
end # function params_intervals

# fitting function
function params_intervals_one_side(
    params::Vector{Float64},
    optim_func::Function,
    constraints_func::Function,
    minmax::Symbol;

    fit_alg::Symbol = :LN_AUGLAG,
    local_alg::Symbol = :LN_NELDERMEAD,
    bounds::Vector{Tuple{Float64,Float64}} = fill((-Inf, Inf), length(params)),
    max_iter::Int64 = 100000,
    ftol_loc::Float64 = 1e-3,
    # ftol_glob::Float64 = 1e-3, # tolerance of global method
    # tol_const::Float64 = 1e-3 # tolerance of constraints
)
    # dim of the problem
    n_params = length(params)

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
    (optf, optx, ret) = optimize(opt, params)
end # function params_intervals_one_side
