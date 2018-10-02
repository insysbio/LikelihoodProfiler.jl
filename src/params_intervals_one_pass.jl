# Calculate interval with ONE_PASS method
function interval_calc(input::ParamInput, ::Val{:ONE_PASS})

    @unpack init_params, id, loss_crit, loss_func, logscale, scan_bound,
    bounds, local_alg, max_iter, ptol, losstol = input

    # set counter scope
    counter::Int64 = 0

    # Output
    result = ParamInterval(
        input,
        Vector{Float64}(2),
        Vector{Symbol}(2),
        Vector{Int64}(2),
        Vector{Float64}(2),

        []
    )

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
            #id,
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
            # result.interval[int_id] = ungarmonize(optf, logscale[id])
            result.ret_codes[int_id] = :BOUNDS_REACHED
        else
            result.interval[int_id] = ungarmonize(optf, logscale[id])
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
    #id::Int64,
    min_max::Symbol;

    fit_alg::Symbol = :LN_AUGLAG,
    local_alg::Symbol = :LN_NELDERMEAD,
    bounds::Vector{Vector{Float64}} = fill(
        [-Inf, Inf], length(init_params)
    ),
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
        #=
        xtol = fill(0., n_params)
        xtol[id] = ftol_loc
        xtol_abs!(local_opt,xtol)
        =#
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
