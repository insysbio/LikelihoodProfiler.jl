# Julia pkg to check parameters identifiability, calculate intervals and plot profiles
#module ParametersIdentification

using NLopt, Plots

#export params_intervals, params_plot


"""
function params_intervals(init_params::Vector{Float64},
                          id::Int64,
                          maxf::Float64,
                          loss_func::Function;
                          logscale::Vector{Bool} = fill(false, length(init_params)),
                          fit_alg::Symbol=:LN_AUGLAG,
                          local_alg::Symbol=:LN_NELDERMEAD,
                          bounds_params::Vector{Vector{Float64}} = fill([-Inf,Inf], length(init_params)),
                          bounds_id::Vector{Float64} = [1e-9, 1e9],
                          max_iter::Float64 = 1e5,
                          tol_glob::Float64 = 1e-1,
                          tol_loc::Float64 = 1e-3,
                          constraints_type::Symbol = :inequality,
                          solver::Symbol=:NLOPT)

# Input:
        init_params - initial parameters vector
        id - id of the parameter for analysis
        maxf - loss function maximum value, "identifiability level"
        loss_func - loss function
        logscale - bool vector length(init_params) where true - log scale / false - direct scale
        fit_alg - fitting algorithm (default - :LN_AUGLAG)
        local_alg - local fitting algorithm (default - :LN_NELDERMEAD)
        bounds_params - bound constraints for all parameters except id
        bounds_id - search bounds for id parameter (default - [1e-9,1e9])
        tol_glob - fitting tolerance for global optimizer (default - NaN)
        tol_loc - fitting tolerance for local optimizer (default - NaN)
        constraints_type - :inequality or :equality constraints
        solver - fitting solver (default - :NLOPT)

# Return:
        confidence intervals evaluation:
        (interval, termination reason, numer of evaluations, loss value)
"""
function params_intervals(init_params::Vector{Float64},
                          id::Int64,
                          maxf::Float64,
                          loss_func::Function;
                          logscale::Vector{Bool} = fill(false, length(init_params)),
                          fit_alg::Symbol=:LN_AUGLAG,
                          local_alg::Symbol=:LN_NELDERMEAD,
                          bounds_params::Vector{Vector{Float64}} = fill([-Inf,Inf], length(init_params)),
                          bounds_id::Vector{Float64} = [1e-9, 1e9],
                          max_iter::Float64 = 1e5,
                          tol_glob::Float64 = NaN,
                          tol_loc::Float64 = NaN,
                          constraints_type::Symbol = :equality,
                          solver::Symbol=:NLOPT) # currently NLOPT is the only option

        # Iterations count
        global count = 0

        # Output
        intervals = Vector{Float64}(2)
        ret_codes = Vector{Symbol}(2)
        count_evals = Vector{Int64}(2)
        loss_final = Vector{Float64}(2)

        # Checks
        (loss_func(init_params) > maxf) && throw(ArgumentError("Check init_params and maxf: loss_func(init_params) should be <= maxf"))
        (init_params[id]<=minimum(bounds_id) || init_params[id]>=maximum(bounds_id)) && throw(ArgumentError("init values are outside of the bounds $bounds_id"))

        # Logscale cheks - parameters
        params = logscale_check_log.(init_params, logscale)

        # Logscale cheks - bounds
        bounds_params[id] = bounds_id
        bounds = logscale_check_bounds.(bounds_params, logscale)

        # Objective function
        optim_func(x, g) = logscale_check_exp(x[id], logscale[id])

        # Constraints function
        function constraints_func(x, g)

                log_chked_x = logscale_check_exp.(x, logscale)
                loss = loss_func(log_chked_x) - maxf

                global count
                count::Int64 += 1

                if (x[id]<=minimum(bounds[id]) || x[id]>=maximum(bounds[id])) && loss < 0.
                        throw(ForcedStop())
                else
                        return loss
                end
        end

        # Confidence interval search
        for minmax in (:min,:max)
                int_id = minmax == :min ? 1 : 2

                (minf,minx,ret) = fitting_params(params,
                                                 optim_func,
                                                 constraints_func,
                                                 id,
                                                 minmax;
                                                 bounds=bounds,
                                                 fit_alg=fit_alg,
                                                 local_alg=local_alg,
                                                 tol_glob=tol_glob,
                                                 tol_loc=tol_loc,
                                                 constraints_type=constraints_type,
                                                 max_iter=max_iter)
                # if bounds reached
                if ret == :FORCED_STOP
                        intervals[int_id] = minmax == :min ? minimum(bounds_id) : maximum(bounds_id)
                        ret_codes[int_id] = :BOUNDS_REACHED
                        count_evals[int_id] = count
                        loss_final[int_id] = loss_func(logscale_check_exp.(minx, logscale))
                else
                        intervals[int_id] = logscale_check_exp(minx[id],logscale[id])
                        ret_codes[int_id] = ret
                        count_evals[int_id] = count
                        loss_final[int_id] = loss_func(logscale_check_exp.(minx, logscale))
                end

                println("id=$id, interval[$int_id] = $(intervals[int_id]), ret_codes=$(ret_codes[int_id]), counts=$(count_evals[int_id])")
                global count = 0
        end

        println(intervals, ", ", ret_codes, ", ", count_evals, ", ", loss_final)
        intervals, ret_codes, count_evals, loss_final
end


# fitting function
function fitting_params(params::Vector{Float64},
                        optim_func::Function,
                        constraints_func::Function,
                        id::Int64,
                        minmax::Symbol;
                        bounds::Vector{Tuple{Float64,Float64}}=fill((-Inf,Inf), length(params)),
                        fit_alg::Symbol=:LN_AUGLAG,
                        local_alg::Symbol=:LN_NELDERMEAD,
                        tol_glob::Float64=1e-1,
                        tol_loc::Float64=1e-3,
                        constraints_type::Symbol = :inequality,
                        max_iter::Float64 = 1e5)

        # dim of the problem
        n_params=length(params)

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
        if fit_alg in (:LN_AUGLAG,:LN_AUGLAG_EQ)
            # local optimizer
            local_opt = Opt(local_alg, n_params)
            # tolerances
            ~isnan(tol_loc) && ftol_abs!(opt,tol_glob)
            ~isnan(tol_loc) && ftol_abs!(local_opt,tol_loc)
            # bound constraints
            lower_bounds!(opt,lb)
            upper_bounds!(opt,ub)
            # adding local optimizer
            local_optimizer!(opt, local_opt)
        else  # other Algorithm
            lower_bounds!(opt,lb)
            upper_bounds!(opt,ub)
            ftol_abs!(opt,tol_glob)
        end

        # max function calls
        maxeval!(opt, Int(max_iter))

        # equality/inequality constraints
        if constraints_type == :equality
                equality_constraint!(opt, constraints_func,tol_loc)
                println("equality constraints selected")
        else
                inequality_constraint!(opt, constraints_func,tol_loc)
                println("inequality constraints selected")
        end
#        println("ftolglob=$tol_glob", )
        # return
        (minf,minx,ret) = optimize(opt,params)
end


"""
    !!! NOT UPDATED !!!
    params_plot(init_params::Vector{Float64},
                id::Int64,
                interval::Tuple{Float64,Float64},
                loss_func::Function,
                maxf::Float64;
                fit_alg::Symbol=:LN_NELDERMEAD,
                tol_glob::Float64=1e-4)

# Input:
        init_params - initial parameters vector
        id - id of the parameter for analysis
        interval - interval for plot
        loss_func - loss function
        maxf - loss function maximum value, "identifiability level"
        fit_alg - fitting algorithm (default - :LN_NELDERMEAD)
        tol_glob - fitting tolerance (default - 1e-4)

# Return:
        parameter profile plot
"""
function params_plot(init_params::Vector{Float64},
                     id::Int64,
                     interval::Tuple{Float64,Float64},
                     loss_func::Function,
                     maxf::Float64;
                     fit_alg::Symbol=:LN_NELDERMEAD,
                     tol_glob::Float64=1e-3)
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
                        ftol_abs!(opt,tol_glob)
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

logscale_check_bounds(b::Vector{Float64}, logscale::Bool) = (b[1] == -Inf) ? (-Inf,log10(b[2])) : (logscale ? Tuple(log10.(b)) : Tuple(b))


#end #module
