struct CICOProfiler2{O,K,A} <: AbstractProfilerMethod
    optimizer::O
    optimizer_opts::K
    auglag_opts::A
    resolution::Int
end

function CICOProfiler2(;
    optimizer,
    optimizer_opts = NamedTuple(),
    auglag_opts = NamedTuple(),
    resolution::Int = 1,
)
    opts = merge((ϵ = 1e-4, ρ_init = 0.1), auglag_opts)
    return CICOProfiler2(optimizer, optimizer_opts, opts, resolution)
end

get_optimizer(method::CICOProfiler2) = method.optimizer
get_optimizer_opts(method::CICOProfiler2) = method.optimizer_opts
get_auglag_opts(method::CICOProfiler2) = method.auglag_opts

mutable struct CICOProfiler2Cache{O,StatsType,R} <: AbstractProfilerMethodCache
    opt_cache::O
    stats::StatsType
    retcode_original::R
end

function init_method_cache(plprob::ProfileLikelihoodProblem, target::AbstractProfileTarget, method::CICOProfiler2, idx, dir, args...; kwargs...)

    optprob_auglag = build_optprob_auglag(plprob, idx, dir)
    opt_cache = SciMLBase.init(optprob_auglag, AugLag(inner = get_optimizer(method); get_auglag_opts(method)...); get_optimizer_opts(method)...)
    return CICOProfiler2Cache(opt_cache, SciMLBase.OptimizationStats(), SciMLBase.ReturnCode.Default)
end

function profiler_step!(branch_cache::ProfileBranchCache, target::AbstractProfileTarget, method_cache::CICOProfiler2Cache)
    @unpack opt_cache = method_cache
    sol = solve_opt_cache(opt_cache)
    x = evaluate_target_f(branch_cache, sol.u)

    if SciMLBase.successful_retcode(sol.retcode)
        branch_cache.pars_cur .= sol.u
        branch_cache.x_cur = x
        branch_cache.obj_cur = evaluate_obj(branch_cache, sol.u)
        branch_cache.iter += 1
    else
        @warn "Solver returned $(sol.retcode) retcode at proposed profile point x = $x. Profiling is interrupted."
    end
    method_cache.retcode_original = sol.retcode
    return nothing
end

function reached_threshold(branch_cache::ProfileBranchCache, method_cache::CICOProfiler2Cache)
   @unpack obj_cur = branch_cache
   obj_level = get_obj_level(branch_cache)
   return abs(obj_cur - obj_level) < method_cache.opt_cache.opt.ϵ_primal
end

function reached_profile_bound(pc::ProfileBranchCache, mc::CICOProfiler2Cache)
    @unpack x_cur, profile_range, obj_cur = pc
    obj_level = get_obj_level(pc)
    tol = mc.opt_cache.opt.ϵ_primal
    bound = isleft(pc) ? profile_range[1] : profile_range[2]

    return isfinite(x_cur) && isfinite(obj_cur) &&
           abs(x_cur - bound) ≤ tol &&
           obj_cur ≤ obj_level + tol
end

function estimate_endpoint(branch_cache::ProfileBranchCache, ::CICOProfiler2Cache)
    return float(branch_cache.x_cur)
end

function build_optprob_auglag(plprob::ProfileLikelihoodProblem, idx, dir)

  optprob = plprob.optprob
  optpars = plprob.optpars

  obj0 = evaluate_obj(plprob, optpars)
  threshold = plprob.threshold
  obj_level = obj0 + threshold

  optf_auglag = build_optf_auglag(plprob, idx, dir)
  lb, ub = bounds_auglag(plprob, idx)

  # Constraints for objective and target profile function
  lcons = [-Inf, get_idx_profile_lb(plprob.target, idx)]
  ucons = [obj_level, get_idx_profile_ub(plprob.target, idx)]

  return OptimizationProblem(optf_auglag, copy(optpars), optprob.p;
   lb, ub, lcons, ucons)
end


function build_optf_auglag(plprob::ProfileLikelihoodProblem, idx, dir)
    opt_f = plprob.optprob.f
    obj_f, obj_f_kwargs = build_objf_auglag(plprob.target, idx, dir)
    cons_f, cons_f_kwargs = build_consf_auglag(plprob.target, opt_f, idx)

    return OptimizationFunction(obj_f, opt_f.adtype; cons = cons_f, obj_f_kwargs..., cons_f_kwargs...)
end

function build_objf_auglag(target::ParameterTarget, idx, dir)
    profile_f = get_profile_f(target, idx)
    obj_f = (x, p) -> -dir * profile_f(x, p)
    kwargs = Dict{Symbol,Any}()

    kwargs[:grad] = function (G, x, p)
        fill!(G, zero(eltype(G)))
        G[idx] = -dir
        return nothing
    end

    kwargs[:hess] = function (H, x, p)
        fill!(H, zero(eltype(H)))
        return nothing
    end
    return obj_f, kwargs
end

function build_objf_auglag(target::FunctionTarget, idx, dir)
    profile_f = get_profile_f(target, idx)
    obj_f = (x, p) -> -dir * profile_f(x, p)
    kwargs = Dict{Symbol,Any}()

    for name in (:grad, :hess)
        derivative = build_derivative_auglag(profile_f, name)
        isnothing(derivative) && continue

        kwargs[name] = function (dest, args...)
            derivative(dest, args...)
            dest .*= -dir
            return nothing
        end
    end

    if !isnothing(profile_f.hv)
        kwargs[:hv] = function (Hv, x, v, p)
            if applicable(profile_f.hv, Hv, x, v, p)
                profile_f.hv(Hv, x, v, p)
            else
                copyto!(Hv, profile_f.hv(x, v, p))
            end
            Hv .*= -dir
            return nothing
        end
    end

    if !isnothing(profile_f.hess_prototype)
        kwargs[:hess_prototype] = copy(profile_f.hess_prototype)
    end
    if !isnothing(profile_f.hess_colorvec)
        kwargs[:hess_colorvec] = copy(profile_f.hess_colorvec)
    end

    return obj_f, kwargs
end

function build_consf_auglag(target::AbstractProfileTarget, opt_f::OptimizationFunction, idx)

    profile_f = get_profile_f(target, idx)

    cons_f = function (res, x, p)
        res[1] = opt_f(x, p)
        res[2] = profile_f(x, p)
        return nothing
    end
    kwargs = Dict{Symbol,Any}()

    obj_grad! = build_derivative_auglag(opt_f, :grad)
    target_grad! = build_target_grad(target, idx)

    # Supply only complete callbacks. If either component is missing,
    # Optimization.jl differentiates both constraints using opt_f.adtype.
    # Both functions must therefore support that backend, even if one has grad.
    if !isnothing(obj_grad!) && !isnothing(target_grad!)

        kwargs[:cons_j] = function (J, x, p)
            obj_grad!(view(J, 1, :), x, p)
            target_grad!(view(J, 2, :), x, p)
            return nothing
        end

        kwargs[:cons_jvp] = function (Jv, x, v, p)
            G = similar(x)
            obj_grad!(G, x, p)
            Jv[1] = dot(G, v)
            target_grad!(G, x, p)
            Jv[2] = dot(G, v)
            return nothing
        end

        kwargs[:cons_vjp] = function (Jv, x, v, p)
            G = similar(x)
            obj_grad!(Jv, x, p)
            Jv .*= v[1]
            target_grad!(G, x, p)
            Jv .+= v[2] .* G
            return nothing
        end
    end

    obj_hess! = build_derivative_auglag(opt_f, :hess)
    target_hess! = build_target_hess(target, idx)

    if !isnothing(obj_hess!) && !isnothing(target_hess!)
        kwargs[:cons_h] = function (Hs, x, p)
            obj_hess!(Hs[1], x, p)
            target_hess!(Hs[2], x, p)
            return nothing
        end
    end

    loss_prototype = opt_f.hess_prototype
    target_prototype = target isa FunctionTarget ? profile_f.hess_prototype : nothing
    if !isnothing(loss_prototype) || !isnothing(target_prototype)
        template = isnothing(loss_prototype) ? target_prototype : loss_prototype
        loss_H = isnothing(loss_prototype) ? zeros(eltype(template), size(template)) : copy(loss_prototype)
        target_H = isnothing(target_prototype) ? zeros(eltype(template), size(template)) : copy(target_prototype)
        kwargs[:cons_hess_prototype] = AbstractMatrix[loss_H, target_H]
    end
    # Let the backend derive constraint colors: a missing component may use a
    # dense fallback, whose pattern need not match either supplied color vector.

    return cons_f, kwargs
end

function build_target_grad(target::ParameterTarget, idx)
    return function (G, x, p)
        G .= 0
        G[idx] = 1
        return nothing
    end
end

function build_target_grad(target::FunctionTarget, idx)
    return build_derivative_auglag(get_profile_f(target, idx), :grad)
end

function build_target_hess(target::ParameterTarget, idx)
    return function (H, x, p)
        H .= 0
        return nothing
    end
end

function build_target_hess(target::FunctionTarget, idx)
    return build_derivative_auglag(get_profile_f(target, idx), :hess)
end

# Normalize each scalar derivative independently; OptimizationFunction's iip flag
# does not establish the calling convention of a user-supplied callback.
function build_derivative_auglag(opt_f::OptimizationFunction, name::Symbol)
    derivative = getproperty(opt_f, name)
    isnothing(derivative) && return nothing
    return function (dest, x, p)
        if applicable(derivative, dest, x, p)
            derivative(dest, x, p)
        elseif applicable(derivative, x, p)
            copyto!(dest, derivative(x, p))
        else
            # A Vector/Matrix-only callback may reject a Jacobian row view.
            scratch = name == :grad ? similar(x) : similar(x, length(x), length(x))
            derivative(scratch, x, p)
            copyto!(dest, scratch)
        end
        return nothing
    end
end

bounds_auglag(plprob::ProfileLikelihoodProblem, idx) = bounds_auglag(plprob.target, plprob.optprob, idx)

function bounds_auglag(target::ParameterTarget, optprob, idx)
    lb, ub = copy_bounds_auglag(optprob)
    if !isnothing(lb)
        lb[idx] = -Inf
    end
    if !isnothing(ub)
        ub[idx] = Inf
    end
    return lb, ub
end

bounds_auglag(target::FunctionTarget, optprob, idx) = copy_bounds_auglag(optprob)

function copy_bounds_auglag(optprob)
    if isnothing(optprob.lb) && isnothing(optprob.ub)
        return nothing, nothing
    end
    n = length(optprob.u0)
    # Fill a missing side so AugLag forwards one-sided bounds to the inner solver.
    lb = isnothing(optprob.lb) ? fill(-Inf, n) : float.(collect(optprob.lb))
    ub = isnothing(optprob.ub) ? fill(Inf, n) : float.(collect(optprob.ub))
    return lb, ub
end
