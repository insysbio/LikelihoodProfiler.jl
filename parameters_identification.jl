# Calculates intervals and profiles for all parameters in minx
using NLopt, Plots

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


        # checks
        (loss_func(init_params) > maxf) && throw(ArgumentError("Check init_params and maxf: loss_func(init_params)should be <= maxf"))
        (init_params[id]<=minimum(init_bounds) || init_params[id]>=maximum(init_bounds)) && throw(ArgumentError("init values are outside of the bounds $init_bounds"))

        # logscale
        params = copy(init_params)
        @. params = logscale_check_log(init_params, logscale)
        bounds = logscale_check_bounds(init_bounds, logscale[id])

        # Constraints function
#        constraints_func(x,g) = loss_func(x,g,logscale) - maxf
        function constraints_func(x, g)

                loss = loss_func(logscale_check_exp.(x, logscale)) - maxf

                global count
                count::Int64 += 1
                println("f_$count($x), loss=$loss")
#                count <= 10 && println("f_$count($x)")

                if (x[id]<=minimum(bounds) || x[id]>=maximum(bounds)) && loss < 0.
                        throw(ForcedStop())
                else
                        return loss
                end
        end

        # Objective function
        function optim_func(x, g)
#                x = logscale_check_exp(x, logscale)
#                count <= 10 && println(x)
                return logscale_check_exp.(x, logscale)[id]
        end

        minmax = (:min,:max)

        for int_id in eachindex(minmax) # parallelize?

                (minf,minx,ret) = fitting_params(params,
                                    optim_func,
                                    constraints_func,
                                    id,
                                    minmax[int_id];
                                    bounds=bounds,
                                    fit_alg=fit_alg,
                                    local_alg=local_alg,
                                    tol_val=tol_val)
                if ret == :FORCED_STOP
                        intervals[int_id] = int_id == 1 ? minimum(init_bounds) : maximum(init_bounds)
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

# log10.
function logscale_check_log(x::Vector{Float64}, logscale::Vector{Bool})
        x_new = copy(x)
        for i in eachindex(x)
                logscale[i] && (x_new[i] = log10(x[i]))
        end
        return x_new
end

logscale_check_log(x::Float64, logscale::Bool) = logscale ? log10(x) : x

# exp10.
function logscale_check_exp(x::Vector{Float64}, logscale::Vector{Bool})
        x_new = copy(x)
        for i in eachindex(x)
                logscale[i] && (x_new[i] = exp10(x[i]))
        end
        return x_new
end

logscale_check_exp(x::Float64, logscale::Bool) = logscale ? exp10(x) : x

# get bounds tuple
logscale_check_bounds(b::Vector{Float64}, logscale::Bool) = logscale ? Tuple(log10.(b)) : Tuple(b)
#=
function logscale_check_bounds(b::Vector{Float64}, logscale::Bool)
        logscale && (b .= log10.(b))
        return Tuple(b)
end
=#

#=
function optim_func(x::Vector{Float64}, g::Vector{Float64}, id::Int64, minmax::Symbol, bounds, logscale::Vector{Bool})
#        logscale_check_exp!(x, logscale)
        println("o($x)")
        if minmax == :min
#                if x[id] >= minimum(bounds)
                        return x[id] #(x[id]-minimum(bounds))^2
#                else return 0. #minimum(bounds)
#                end
        elseif minmax == :max
#                if x[id] <= maximum(bounds)
                        return x[id] #-(x[id]-maximum(bounds))^2
#                else return 0. #maximum(bounds)
#                end
        end
end
=#
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
        # bounds

        #=
        lb = fill(-Inf, n_params)
        ub = fill(Inf, n_params)
        lb[id] = minimum(bounds)
        ub[id] = maximum(bounds)
        =#

        lb = fill(minimum(bounds), n_params)
        ub = fill(maximum(bounds), n_params)


        # optimization obj
        opt = Opt(fit_alg, n_params)

        # min or max
        if minmax == :min
            min_objective!(opt, optim_func)

        elseif minmax == :max
            max_objective!(opt, optim_func)
        end

        if fit_alg in (:LN_AUGLAG,:LN_AUGLAG_EQ)
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



# to do!
function paramPlots(minx::Vector{Float64}, idx::Int, intervals::Vector; logscale::Bool = false)


        profile_i = Tuple{Float64,Float64}
        InitVals = [minx[j] for j in 1: length(minx) if j != idx]
        minx_i = minx[idx]
        intervals_i = intervals[idx]

        if logscale == true
                minx = log10.(minx)
        end

        param_plot_i = y -> likelihood_i(y, idx, InitVals; logscale=logscale)
        plot_minmax = (min(minx_i,intervals_i[1]) - minx_i/10., max(minx_i,intervals_i[2]) + minx_i/10.)
        profile_i = adapted_grid(param_plot_i,(minimum(plot_minmax),maximum(plot_minmax)); max_recursions = 7)

        profile_i
end

# Calculates loss for fixed parameter minx_i f(minx_i) -> minf_i
function likelihood_i(minx_i::Float64, idx::Int64, InitParams::Vector{Float64}; logscale::Bool=false)
    prob_i = (prob, uParams) -> getMonteCarloProblemFixedParam_i(prob, uParams, idx, minx_i; logscale=logscale)
    obj_i = build_loss_objective(empty_ode_problem,Rosenbrock23(),ofv_loss,prob_generator=prob_i,
    verbose=false, verbose_opt=true, verbose_steps=50, maxiter=100000,
    num_monte=n_sim, saveat=t_save, save_idxs = [2])
    @time (minf_i,nextInit,ret) = fiting_opt(:LN_NELDERMEAD, InitParams, obj_i; logscale=logscale)
    return minf_i
end

function getMonteCarloProblemFixedParam_i(prob, uParams::Array{Float64,1}, idx::Int64, minx_i::Float64; logscale::Bool=false)
    newParams=[uParams[j] for j in eachindex(uParams)]
    insert!(newParams, idx, minx_i)
    if logscale == true
            newParams = exp10.(newParams)
    end
    MonteCarloProblem(prob, prob_func=getODEProblemGenerator(newParams))
end

function ode_loss(InitParams::Vector, grad::Vector = []; logscale::Bool = false)
   global count
   count::Int += 1
   println("f_$count($InitParams)")
   constrainedProblem = getMonteCarloProblem(empty_ode_problem, InitParams; logscale=logscale)
   constrainedSol = solve(constrainedProblem, Rosenbrock23(), num_monte=n_sim, saveat=t_save, save_idxs = [2])
   ofv_loss(constrainedSol)
end


# fiting using NLopt
function fiting_func(alg::Symbol,
                     init_params::Vector{Float64},
                     obj::Function;
                     bounds::Tuple{Float64,Float64} = (1e-9,1e9),
                     constraints::Bool = false,
                     constraints_func::Function = x->0.0,
                     local_alg::Symbol = :LN_NELDERMEAD,
                     minmax::Symbol = :min,
                     idx::Int64 = 1,
                     tol_val = 1e-3,
                     logscale::Bool=false)


        n_params=length(init_params)
        opt = Opt(alg, n_params)

        if logscale == true
            bounds = log10.(bounds)
        end

        if minmax == :min
            min_objective!(opt, obj)
            stop_val = maximum(bounds)
#            if constraints == true
#                inequality_constraint!(opt, (x,g)->(-9.0-x[idx]), 1e-5)
#            end
        elseif minmax == :max
            max_objective!(opt, obj)
            stop_val = minimum(bounds)
#            if constraints == true
#                inequality_constraint!(opt, (x,g)->(-9.0+x[idx]), 1e-5)
#            end
        end


        if alg == :LN_AUGLAG
            local_opt = Opt(local_alg, n_params)
            local_optimizer!(opt, local_opt)
    #        lb = fill(-Inf, n_params)
    #        ub = fill(Inf, n_params)
    #        lb[idx] = minimum(bounds)
    #        ub[idx] = maximum(bounds)
    #        lower_bounds!(opt,lb)
    #        upper_bounds!(opt,ub)
            stopval!(opt, stop_val)
    #        x_tol = fill(0.0, n_params)
    #        x_tol[idx] = 1e-3
    #        xtol_abs!(local_opt,1e-3)
    #        initial_step!(opt, fill(0.2, n_params))
    #        ftol_abs!(local_opt,1e-3)
        else
            lower_bounds!(opt,fill(minimum(bounds), n_params))
            upper_bounds!(opt,fill(maximum(bounds), n_params))
        end

#        xtol_abs!(opt,1e-3)
#        xtol_abs!(opt,1e-3)
        ftol_abs!(opt,tol_val)
        maxeval!(opt, Int(1e5))
        constraints && inequality_constraint!(opt, constraints_func, 1e-5)

        (minf,minx,ret) = optimize(opt,init_params)
        if logscale == true
            minx = exp10.(minx)
        end
    #    println("$(initial_step(opt, uParamsInit))")
        (minf,minx,ret)
end
