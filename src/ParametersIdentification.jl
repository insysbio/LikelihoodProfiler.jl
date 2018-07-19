# Julia pkg to check parameters identifiability, calculate intervals and plot profiles
module ParametersIdentification

using NLopt, Plots

export params_intervals, params_plot


"""
    params_intervals(init_params::Vector{Float64},
            id::Int64,
            maxf::Float64,
            loss_func::Function,
            logscale::Vector{Bool};
            fit_alg::Symbol=:LN_AUGLAG,
            local_alg::Symbol=:LN_SBPLX,
            init_bounds::Vector{Float64} = [1e-9,1e9],
            tol_val::Float64 = 1e-4,
            solver::Symbol=:NLOPT)

# Input:
        init_params - initial parameters vector
        id - id of the parameter for analysis
        maxf - loss function maximum value, "identifiability level"
        loss_func - loss function
        logscale - bool vector length(init_params) where true - log scale / false - direct scale
        fit_alg - fitting algorithm (default - :LN_AUGLAG)
        local_alg - local fitting algorithm (default - :LN_SBPLX)
        init_bounds - search bounds (default - [1e-9,1e9])
        tol_val - fitting tolerance (default - 1e-4)
        solver - fitting solver (default - :NLOPT)

# Return:
        confidence intervals evaluation:
        (interval, termination reason, numer of evaluations)
"""
function params_intervals(init_params::Vector{Float64},
                          id::Int64,
                          maxf::Float64,
                          loss_func::Function,
                          logscale::Vector{Bool};
                          fit_alg::Symbol=:LN_AUGLAG,
                          local_alg::Symbol=:LN_SBPLX,
                          init_bounds::Vector{Float64} = [1e-9,1e9],
                          tol_val::Float64 = 1e-4,
                          solver::Symbol=:NLOPT)

        # Iterations count
        global count = 0

        # Output
        intervals = Vector{Float64}(2)
        ret_codes = Vector{Symbol}(2)
        count_evals = Vector{Int64}(2)


        # Checks
        (loss_func(init_params) > maxf) && throw(ArgumentError("Check init_params and maxf: loss_func(init_params)should be <= maxf"))
        (init_params[id]<=minimum(init_bounds) || init_params[id]>=maximum(init_bounds)) && throw(ArgumentError("init values are outside of the bounds $init_bounds"))

        # Logscale
        params = copy(init_params)
        @. params = logscale_check_log(init_params, logscale)
        bounds = logscale_check_bounds(init_bounds, logscale[id])

        # Objective function
        optim_func(x, g) = logscale_check_exp(x[id], logscale[id])

        # Constraints function
        function constraints_func(x, g)

                loss = loss_func(logscale_check_exp.(x, logscale)) - maxf

                global count
                count::Int64 += 1
                println("f_$count($x), loss=$loss")

                if (x[id]<=minimum(bounds) || x[id]>=maximum(bounds)) && loss < 0.
                        throw(ForcedStop())
                else
                        return loss
                end
        end


        for minmax in (:min,:max) # parallelize?
                int_id = minmax == :min ? 1 : 2

                (minf,minx,ret) = fitting_params(params,
                                    optim_func,
                                    constraints_func,
                                    id,
                                    minmax;
                                    bounds=bounds,
                                    fit_alg=fit_alg,
                                    local_alg=local_alg,
                                    tol_val=tol_val)
                if ret == :FORCED_STOP
                        intervals[int_id] = minmax == :min ? minimum(init_bounds) : maximum(init_bounds)
                        ret_codes[int_id] = :BOUNDS_REACHED
                        count_evals[int_id] = count
                else
                        intervals[int_id] = logscale_check_exp(minx[id],logscale[id])
                        ret_codes[int_id] = ret
                        count_evals[int_id] = count
                end

                global count = 0
                println("interval[$int_id] = $(intervals[int_id])")
        end

        intervals, ret_codes, count_evals
end


# fitting function

function fitting_params(params::Vector{Float64},
                        optim_func::Function,
                        constraints_func::Function,
                        id::Int64,
                        minmax::Symbol;
                        bounds::Tuple{Float64,Float64}=(1e-9,1e9),
                        fit_alg::Symbol=:LN_AUGLAG,
                        local_alg::Symbol=:LN_NELDERMEAD,
                        tol_val::Float64=1e-3)

        # dimention of the problem
        n_params=length(params)

        println("init params == $params")
        println("bounds == $bounds")
        println("minmax == $minmax")


        # optimization obj
        opt = Opt(fit_alg, n_params)

        # min or max
        if minmax == :min
            min_objective!(opt, optim_func)

        elseif minmax == :max
            max_objective!(opt, optim_func)
        end

        if fit_alg in (:LN_AUGLAG,:LN_AUGLAG_EQ)
            #=
            lb = fill(-Inf, n_params)
            ub = fill(Inf, n_params)
            lb[id] = minimum(bounds)
            ub[id] = maximum(bounds)
            =#

            lb = fill(minimum(bounds), n_params)
            ub = fill(maximum(bounds), n_params)
            local_opt = Opt(local_alg, n_params)
            ftol_abs!(opt,tol_val)
#            lower_bounds!(opt,lb)
#            upper_bounds!(opt,ub)
            local_optimizer!(opt, local_opt)
        #    stopval!(local_opt, stop_val)

        else
            lower_bounds!(opt,fill(minimum(bounds), n_params))
            upper_bounds!(opt,fill(maximum(bounds), n_params))
        end

#        println(opt.cb[1])

#        ftol_abs!(opt,tol_val)

        maxeval!(opt, Int(1e5))
        inequality_constraint!(opt, constraints_func,tol_val)
#        inequality_constraint!(opt, (x,g)->(x[id]-1e9))
#        inequality_constraint!(opt, (x,g)->(1e9-x[id]))
        (minf,minx,ret) = optimize(opt,params)
end

"""
    params_plot(init_params::Vector{Float64},
                id::Int64,
                interval::Tuple{Float64,Float64},
                loss_func::Function,
                maxf::Float64;
                fit_alg::Symbol=:LN_NELDERMEAD,
                tol_val::Float64=1e-4)

# Input:
        init_params - initial parameters vector
        id - id of the parameter for analysis
        interval - interval for plot
        loss_func - loss function
        maxf - loss function maximum value, "identifiability level"
        fit_alg - fitting algorithm (default - :LN_NELDERMEAD)
        tol_val - fitting tolerance (default - 1e-4)

# Return:
        parameter profile plot
"""
function params_plot(init_params::Vector{Float64},
                     id::Int64,
                     interval::Tuple{Float64,Float64},
                     loss_func::Function,
                     maxf::Float64;
                     fit_alg::Symbol=:LN_NELDERMEAD,
                     tol_val::Float64=1e-4)
        # cheks
        (loss_func(init_params) > maxf) && throw(ArgumentError("Check init_params and maxf: loss_func(init_params)should be <= maxf"))

        params = copy(init_params)

        function profile_func(x::Float64)
                if length(params) == 1
                        return loss_func([x])
                else
                        fit_params_func = (p,g) -> loss_func(insert!(copy(p),id,x))
                        opt = Opt(fit_alg, length(init_params)-1)
                        min_objective!(opt, fit_params_func)
                        ftol_abs!(opt,tol_val)
                        (loss,minx,ret) = optimize(opt,deleteat!(copy(params),id))
                        return loss
                end
        end
        grid = adapted_grid(profile_func,interval; max_recursions = 7)

        gr = plot(grid[1], interval, grid[2], label="param1");
        hline!([maxf], label="alpha_level $maxf")
        scatter!(grid[1], grid[2], label="plot_grid")
#        scatter!(intervals[1],[alpha,alpha], label="alpha_level intersections")

        return gr
end

# Logscale checks
logscale_check_log(x::Float64, logscale::Bool) = logscale ? log10(x) : x

logscale_check_exp(x::Float64, logscale::Bool) = logscale ? exp10(x) : x

logscale_check_bounds(b::Vector{Float64}, logscale::Bool) = logscale ? Tuple(log10.(b)) : Tuple(b)


end #module
