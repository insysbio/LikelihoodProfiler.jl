# Pkg.add("NLopt")

using NLopt

function params_intervals_d2d(
    init_params::Vector{Float64},
    id::Int64,
    loss_crit::Float64,
    loss_func::Function;

    logscale_all::Bool = false,
    logscale::Vector{Bool} = fill(logscale_all, length(init_params)),
    scan_bound::Vector{Float64} = ungarmonize.(
        [-9., 9.],
        logscale[id]
    ),
    fit_alg::Symbol = :LN_BOBYQA, # not used
    local_alg::Symbol = :LN_NELDERMEAD,
    bounds::Vector{Vector{Float64}} = ungarmonize.( # not used
        fill([-Inf, Inf], length(init_params)),
        logscale
    ),
    max_iter::Int64 = 100000, # not used
    ftol_loc::Float64 = 1e-3,

    q::Float64 = ftol_loc # d2d parameter
)
    # set counter scope
    counter::Int64 = 0

    # Output
    intervals = Vector{Float64}(2)
    ret_codes = Vector{Symbol}(2)
    count_evals = Vector{Int64}(2)
    loss_final = Vector{Float64}(2)

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

    function loss_func_upd(params::Vector{Float64})
        loss = loss_func(params)
        counter += 1

        # println("f_$counter($params), loss=$loss")
        return loss
    end

    # Logscale
    params = garmonize.(init_params, logscale)
    scan_bound_garm = garmonize(scan_bound, logscale[id])

    function profile_func(id::Int64, params::Vector{Float64})
        if length(params) == 1
            loss = loss_func_upd(params)
            return params, loss
        else
            function fit_params_func(p,g)
                loss = loss_func_upd(insert!(copy(p), id, params[id]))
                loss
             end
            opt = Opt(local_alg, length(params)-1)
            min_objective!(opt, fit_params_func)
            ftol_abs!(opt, ftol_loc)
            (loss, minx, ret) = optimize(opt, deleteat!(copy(params), id))
            return insert!(minx,id,params[id]), loss
        end
    end

    # init d2d settings
    init_loss = loss_func_upd(params)
    delta_alpha = loss_crit - init_loss
    q_delta_alpha = q * delta_alpha
    dps = zeros(length(params))
    minstepsize = 1e-6 # minumum size of a step
    maxstepsize = 0.4 * params[id] # XXX: magic number
    init_step = maxstepsize / 2.1 # XXX: magic number
    # maxsteps = 1e3 # not implemented

    for step in (-init_step, init_step)
        # Iterations counter
        counter = 0

        step < 0. ? p_id = 1 : p_id = 2
        ps = copy(params)
        dps[id] = step
        loss = init_loss
        while true
            step = getStepDirect(
                id,
                ps,
                step,
                q_delta_alpha,
                loss_func_upd,
                loss,
                scan_bound_garm,
                minstepsize,
                maxstepsize
            )
            if isnan(step)
                    ret_codes[p_id] = :BOUNDS_REACHED
                    break
            end
            dps[id] = step
            p_trial = ps + dps
            ps, loss = profile_func(id, p_trial)
            if loss > init_loss + delta_alpha #*1.2
                    ret_codes[p_id] = :FTOL_REACHED
                    break
            end
        end
        intervals[p_id] = ps[id]
        count_evals[p_id] = counter
        loss_final[p_id] = loss
    end
    println("PLE intervals = $intervals")

    return intervals, ret_codes, count_evals, loss_final
end # function

function getStepDirect(
    id::Int64,                # p[jk] id
    ps::Vector{Float64},      # current p values
    last_step::Float64,       # current dp
    q_delta_alpha::Float64,
    loss_func_upd::Function,
    loss::Float64,            # f threshold
    bounds,
    minstepsize::Float64,
    maxstepsize::Float64
)
    # bounds
    lb = bounds[1]
    ub = bounds[2]

    # steps
    #    minstepsize = 1e-6 # minumum size of a step
    #    maxstepsize = 0.2*ps[id] # maximum size of a step
    stepfaktor = 2

    step = last_step

    if (ps[id]+step <= lb + minstepsize || ps[id]+step >= ub - minstepsize) # jk hit boundaries
        step = step / stepfaktor
        if abs(step) < minstepsize
            if (step>0)
                lbub = "upper"
            else
                lbub = "lower"
            end
        end
        println("PLE_$id parameter hit boundary")
        return NaN
    else
        ps[id] = ps[id] + step
        loss_trial = loss_func_upd(ps)
        if (loss_trial - loss) > q_delta_alpha
            step = step / stepfaktor
            if abs(step) < minstepsize
                println("WARNING: could not control step size")
                step = step * stepfaktor
            end
        else
            if ((loss_trial - loss) < q_delta_alpha) && ((loss_trial - loss) >= 0.)
                step = step * stepfaktor
                if abs(step) > maxstepsize
                    step = maxstepsize*sign(step)
                end
            end
        end

        println("step=$step")
        return step
    end
end
