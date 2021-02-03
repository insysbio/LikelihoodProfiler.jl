### !!! Current version of auglag code

using NLopt, LinearAlgebra

# structure representing NLModel
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

function plmodel(loss_func, x0; kwargs...)
  [plmodel(i, loss_func,x0; kwargs...) for i in eachindex(x0)]
end

function plmodel(id::Int, loss_func, x0;
  lbounds=fill(-1e9,length(x0)), ubounds=fill(1e9,length(x0)), scan_range=(-100,100), crit_level=9.0)

  @assert loss_func(x0)<crit_level "Check x0 and crit_level: loss_func(x0) should be < crit_level"
  @assert lbounds[id] < scan_range[1] && scan_range[2] < ubounds[id] "scan_range should be inside of feasible are. Check lbounds and ubounds for $id parameter"
  
  x0 .= promote_type(x0)
  T = eltype(x0)

  PLModel(
    x->x[id],
    loss_func,
    x0,
    eltype(lbounds) == T ? lbounds : T.(lbounds),
    eltype(ubounds) == T ? ubounds : T.(ubounds),
    T.(scan_range),
    typeof(crit_level) ==T ? crit_level : T(crit_level)
  )
end

# structure holding Lagrangian related data
mutable struct Auglag{F1,F2,T}
    scan_func::F1
    constr_func::F2
    #loss_func::F2
    #crit_level::T
    λ::T
    μ::T
    #direction::Int
end

obj(auglag::Auglag, direction) = (x,g) -> direction*auglag.scan_func(x) - 
  auglag.λ*auglag.constr_func(x) + 
    0.5*auglag.μ*(auglag.constr_func(x))^2

obj_ineq(auglag::Auglag, direction) = (x,g) -> direction*auglag.scan_func(x) - 
    auglag.λ*auglag.constr_func(x) + 
      0.5*auglag.μ*(auglag.constr_func(x))^2

init_λ!(auglag::Auglag, λ) = auglag.λ = λ
init_μ!(auglag::Auglag, μ) = auglag.μ = μ
update_λ!(auglag::Auglag, x) = auglag.λ = auglag.λ-auglag.μ*auglag.constr_func(x)
update_μ!(auglag::Auglag) = auglag.μ = 100*auglag.μ


identify_auglag(model_vec::Vector{M}; kwargs...) where M <: PLModel = 
  [identify_auglag(m; kwargs...) for m in model_vec]

function identify_auglag(
  model::PLModel;
  subsolver::Symbol = :LN_NELDERMEAD,
  max_time::Real = 30.0,
  max_eval::Int = 10000,
  #atol::Real = 1e-8, 
  #rtol::Real = 1e-8, 
  converge_tol::Real = 1e-3,
  constr_tol::Real = 1e-3
)   

  lbounds = model.lbounds
  ubounds = model.ubounds

  x0 = model.x0
  x0 .= max.(lbounds, min.(x0, ubounds))

  # Final tolerances
  η_final = constr_tol
  ω_final = converge_tol

  # Starting values as in Nocedal p.520
  μ₀ = 10.
  ω₀ = ω_final # 1/μ₀ # not implemented
  η₀ = 1/μ₀^0.1

  # Init lagrange multiplier
  λ₀ = 1.0

  # Init Augmented Lagrangian
  constr_func(x) = model.loss_func(x) - model.crit_level
  auglag = Auglag(model.scan_func, constr_func, λ₀, μ₀)

  ep_lower = identify_endpoint(auglag, copy(x0), :lower, η₀, ω₀, η_final, ω_final, model.scan_range;
    alg=subsolver, lb=lbounds, ub=ubounds, max_time=max_time, max_eval=max_eval) # can be changed in Julia 1.5.3

  init_λ!(auglag, λ₀)
  init_μ!(auglag, μ₀)

  ep_upper = identify_endpoint(auglag, copy(x0), :upper, η₀, ω₀, η_final, ω_final, model.scan_range;
    alg=subsolver, lb=lbounds, ub=ubounds, max_time=max_time, max_eval=max_eval) # can be changed in Julia 1.5.3

  return (ep_lower, ep_upper)
end

function identify_endpoint(auglag, x, ep_type, η, ω, η_final, ω_final, scan_range; 
  max_time=30.0, max_eval=10000, kwargs...)
  
  # auglag functions
  scan_func = auglag.scan_func
  constr_func = auglag.constr_func

  # search direction
  direction_sign = ep_type == :lower ? 1 : -1

  iter = 0
  evals_tot = 0
  start_time = time()
  el_time = 0.0

  ep = scan_func(x)
  constr_viol = constr_func(x)
  status = :none

  while status ∉ (:identifiable,  :unidentifiable,  :failure,  :tired)
            
    # solve subproblem
    (optf, optx, ret_code, evals) = subproblem_opt(obj(auglag,direction_sign), x; 
      tol=ω, kwargs...)

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

      # Test unidentifiability
      if cf ≤ 0.0 && !(scan_range[1] ≤ scan_val ≤ scan_range[2])
                    
        ep = nothing
        constr_viol = constr_func(x)
        status = :unidentifiable
        break
      end

      # Test identifiability
      if normcf ≤ η && (scan_range[1] ≤ scan_val ≤ scan_range[2])
        # test for convergence
        if normcf ≤ η_final  # && normgp ≤ ϵd gradient projection method not implemented
          ep = scan_val
          constr_viol = constr_func(x)
          status = :identifiable
          break
        end
        # update multipliers, tighten tolerances
        update_λ!(auglag, x)
        η = max(η/auglag.μ^0.9, η_final)
        # ω = ω/auglag.μ # not implemented
      else
        # increase penalty parameter, tighten tolerances
        update_μ!(auglag)
        η = max(η/auglag.μ^0.1, η_final)
        # ω = 1/auglag.μ # not implemented
      end

    else
      ep = nothing
      constr_viol = nothing
      status = opt_status
    end
  end

  return (ep,constr_viol,status,evals_tot,el_time)
end

function subproblem_opt(lagrange_obj, x; 
  lb=fill(1e-9,length(x)), ub=fill(1e9,length(x)), alg=:LN_NELDERMEAD, tol=1e-4, max_eval=50000, max_time=30.0)
    
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

