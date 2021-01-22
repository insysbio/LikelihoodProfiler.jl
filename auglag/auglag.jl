### !!! Current version is  auglag2.jl

using NLopt, LinearAlgebra

struct PLModel{SF,LF,T}
    scan_func::SF
    loss_func::LF
    x0::Vector{T}
    lbounds::Vector{T}
    ubounds::Vector{T}
    scan_range::Tuple{T,T}
    crit_level::T
    #nvar::Int
end

mutable struct Auglag{F1,F2,T}
    scan_func::F1
    constr_func::F2
    λ::T
    μ::T
    #η::T
    #ω::T 
    #η_final::T
    #ω_final::T
    direction::Int
end

update_direction!(auglag::Auglag, direction::Int) = auglag.direction = direction
update_λ!(auglag::Auglag, λ) = auglag.λ = λ
update_μ!(auglag::Auglag, μ) = auglag.μ = μ

obj(auglag::Auglag) = (x,g)-> auglag.direction*auglag.scan_func(x) + auglag.λ*auglag.constr_func(x) + 0.5*auglag.μ*(auglag.constr_func(x))^2


"""
Implementation of an augmented Lagrangian method. The following keyword parameters can be passed:
- μ: Starting value of the penalty parameter (default: 10.0)
- atol: Absolute tolerance used in dual feasibility measure (default: 1e-8)
- rtol: Relative tolerance used in dual feasibility measure (default: 1e-8)
- ctol: (Absolute) tolerance used in primal feasibility measure (default: 1e-8)
- max_iter: Maximum number of iterations (default: 1000)
- max_time: Maximum elapsed time in seconds (default: 30.0)
- max_eval: Maximum number of objective function evaluations (default: 100000)
- subsolver_logger: Logger passed to `tron` (default: NullLogger)
"""
function identify_auglag(
    model::PLModel,
    subsolver::Symbol = :LN_NELDERMEAD;
    μ::Real = 10.0, # penalty parameter
    max_time::Real = 30.0,
    max_eval::Int = 10000,
    #atol::Real = 1e-8, 
    #rtol::Real = 1e-8, 
    converge_tol::Real = 1e-3,
    constr_tol::Real = 1e-3
)   
    T = eltype(model.x0)
    lb = model.lbounds
    ub = model.ubounds

    x0 = model.x0
    x0 .= max.(lb, min.(x0, ub))
    nvar = length(x0)

    scan_func = model.scan_func
    loss_func = model.loss_func
    constr_func(x) = loss_func(x) - model.crit_level

    # Final tolerances
    η_final = constr_tol
    ω_final = converge_tol

    # Starting values as in Nocedal
    μ₀ = 10.
    ω₀ = ω_final #1/μ₀ # not implemented
    η₀ = 1/μ₀^0.1

    # Init lagrange multiplier
    λ₀ = 1.0

    # Init Augmented Lagrangian
    auglag = Auglag(scan_func, constr_func, λ₀, μ₀, 1)

    results = Vector{Any}(undef, 2)
    for (i, endpoint) in enumerate((:lower, :upper))
        
        x = copy(x0)
        μ = μ₀; ω = ω₀; η = η₀

        # Lagrange setup
        direction = sign_of_direction(Val{endpoint})
        λ = λ₀
        update_direction!(auglag, direction)
        update_λ!(auglag, λ)
        update_μ!(auglag, μ)

        iter = 0
        evals_tot = 0
        start_time = time()
        el_time = 0.0

        ep = scan_func(x)
        loss = loss_func(x)
        status = :none

        while status ∉ (:identifiable,  :unidentifiable,  :failure,  :tired)
            
            # solve subproblem
            (optf, optx, ret_code, evals) = subproblem_opt(obj(auglag), x; lb=lb, ub=ub, alg=subsolver, tol=ω, max_eval=max_eval)
            
            if ret_code in (:FTOL_REACHED, :XTOL_REACHED, :SUCCESS)
                opt_status = :success
            elseif ret_code in (:MAXEVAL_REACHED, :MAXTIME_REACHED)
                opt_status = :tired
            else
                opt_status = :failure
            end

            # evals and time measure
            iter += 1
            evals_tot += evals
            el_time = time() - start_time
            opt_status = evals_tot > max_eval || el_time > max_time ? :tired : opt_status

            if opt_status == :success
                # new optimal x
                x .= optx

                # endpoint estimation
                scan_val = scan_func(x)

                # Norm of constr function
                cf = constr_func(x)
                normcf = norm(cf)
               # @show x
               # @show cf
               # @show scan_val
               # @show cf ≤ 0.0
               # @show !(model.scan_range[1] ≤ scan_val ≤ model.scan_range[2])
                # Test unidentifiability
                if cf ≤ 0.0 && !(model.scan_range[1] ≤ scan_val ≤ model.scan_range[2])
                    
                    ep = nothing
                    loss = loss_func(x)
                    status = :unidentifiable
                    break
                end

                # Test identifiability
                if normcf ≤ η
                    # test for convergence
                    if normcf ≤ η_final  # && normgp ≤ ϵd gradient projection method not implemented
                        ep = scan_val
                        loss = loss_func(x)
                        status = :identifiable
                        break
                    end
                    # update multipliers, tighten tolerances
                    update_λ!(auglag, auglag.λ-auglag.μ*constr_func(x))
                    η = max(η/auglag.μ^0.9, η_final)
                    # ω = ω/auglag.μ # not implemented
                else
                    # increase penalty parameter, tighten tolerances
                    update_μ!(auglag,100*auglag.μ)
                    η = max(η/auglag.μ^0.1, η_final)
                    # ω = 1/auglag.μ # not implemented
                end
            else
                ep = nothing
                loss = nothing
                status = opt_status
            end
        end
        
        results[i] = (ep,loss,status,evals_tot,el_time)
    end

    return results
end


function identify_endpoint(auglag::Auglag, x0, )
    
end


function subproblem_opt(lagrange_obj, x; lb=fill(1e-9,length(x)), ub=fill(1e9,length(x)), alg=:LN_NELDERMEAD, tol=1e-4, max_eval=50000)
    
    # NLopt optimization problem 
    nparams = length(x)
    opt = Opt(alg, nparams)
    
    # Box constraints
    opt.lower_bounds = lb
    opt.upper_bounds = ub

    # Tolerance setup - ??
    ftol_abs!(opt, tol)

    maxeval!(opt, max_eval)

    min_objective!(opt, lagrange_obj)

    (optf, optx, ret) = optimize(opt, x)

    return (optf, optx, ret, opt.numevals)
end

sign_of_direction(::Type{Val{:lower}}) = 1
sign_of_direction(::Type{Val{:upper}}) = -1


# test
model_1p = PLModel((x)->x[1], (x)->5.0 + (x[1]-3.0)^2, [3.],[-1e-9],[1e9], (-1e-11,100.), 9.0)
model_2p = PLModel((x)->x[2], (x)->5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2, [3.,4.],[-1e-9,-1e-9],[1e9,1e9], (-1e-11,100.), 9.0)
model_2p = PLModel((x)->x[2], (x)->5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2, [3.,4.],[-1e-9,-1e-9],[1e9,1e9], (-1e-11,100.), 9.0)
model_2p_1im = PLModel((x)->x[2], (x)->5.0 + (x[1]-3.0)^2 + 0*x[2], [3.,4.],[-1e-9,-1e-9],[1e9,1e9], (-1e-11,100.), 9.0)
model_3p_1im = PLModel((x)->x[3], (x)->5.0 + (x[1]-3.0)^2 + (x[2]/x[3]-4.0)^2, [3.,4.,3.0],[-1e-9,-1e-9,1e-9],[1e9,1e9,1e9], (-1e-11,100.), 9.0)