# Calculate interval with D2D_PLE method
function interval_calc(input::ParamInput, ::Val{:D2D_PLE})
    # unpack parameters e.x input.init_params to init_params
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

    function loss_func_upd(params::Vector{Float64})
        loss = loss_func(params)
        counter += 1

        return loss
    end

    # Logscale
    params = garmonize.(init_params, logscale)
    scan_bound_garm = garmonize(scan_bound, logscale[id])


    lb = minimum.(bounds)
    ub = maximum.(bounds)

    function profile_func(id::Int64, params::Vector{Float64})
        if length(params) == 1
            loss = loss_func_upd(params)
            return params, loss
        else
            fit_params_func(p,g) = loss_func_upd(p)

            opt = Opt(local_alg, length(params))
            min_objective!(opt, fit_params_func)
            ftol_abs!(opt, losstol)

            # exclude params[id] from optimization
            lb[id] = params[id]
            ub[id] = params[id]
            lower_bounds!(opt, lb)
            upper_bounds!(opt, ub)

            (loss,minx,ret) = optimize(opt, params)
            return minx, loss
        end
    end

    # init d2d settings
    q = losstol # ?? losstol
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
                    result.ret_codes[p_id] = :BOUNDS_REACHED
                    break
            end
            dps[id] = step
            p_trial = ps + dps
            ps, loss = profile_func(id, p_trial)
            if loss > init_loss + delta_alpha #*1.2
                    result.ret_codes[p_id] = :FTOL_REACHED
                    break
            end
        end
        result.intervals[p_id] = ps[id]
        result.count_evals[p_id] = counter
        result.loss_final[p_id] = loss
    end
    #println("PLE intervals = $intervals")

    result
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
        println("WARNING: PLE_$id parameter hit boundary")
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

        #println("step=$step")
        return step
    end
end
