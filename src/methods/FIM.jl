############################## FIMProfiler ##############################

"""
    FIMProfiler

Fisher Information Matrix (FIM)-based asymptotic confidence intervals (Wald approximation).
By default this method reuses Hessian logic from `OptimizationProblem` (user-supplied Hessian or AD backend).
The confidence interval is computed as `θ̂ ± z * sqrt(Σ[idx, idx])`, where 
        - `θ̂` is the `optpars[idx]`, 
        - `z` is the quantile of the chi-squared distribution corresponding to the `conf_level` and `df` parameters of the `ProfileLikelihoodProblem`,
        - `Σ` is the covariance matrix obtained by inverting the FIM.

`covariance_factor` controls Hessian/objective scaling conventions. Common choices are:
- `1.0` when Hessian is for `-logL`.
- `2.0` when Hessian is for `-2logL` and you want covariance on the `-logL` scale.

Any strictly positive value is allowed (not only `1` or `2`), which can be useful for calibrated or robust variance scaling.

### Fields

- `inversion::Symbol`: Matrix inversion strategy (`:cholesky`, `:pinv`).
- `clamp_to_bounds::Bool`: Clip estimated interval endpoints to profile bounds.
- `covariance_factor::Real`: Multiplicative factor applied to `inv(H)` to obtain covariance (`Σ = covariance_factor * inv(H)`).
"""
Base.@kwdef struct FIMProfiler <: AbstractProfilerMethod
  inversion::Symbol = :cholesky
  clamp_to_bounds::Bool = true
  covariance_factor::Float64 = 1.0

  function FIMProfiler(inversion::Symbol, clamp_to_bounds::Bool, covariance_factor::Real)
    inversion in (:cholesky, :pinv) ||
      throw(ArgumentError("`inversion` must be one of :cholesky, :pinv (got $inversion)."))

    covariance_factor > 0 || throw(ArgumentError("`covariance_factor` must be strictly positive (got $covariance_factor)."))

    new(inversion, clamp_to_bounds, float(covariance_factor))
  end
end

get_inversion(fp::FIMProfiler) = fp.inversion
get_clamp_to_bounds(fp::FIMProfiler) = fp.clamp_to_bounds
get_covariance_factor(fp::FIMProfiler) = fp.covariance_factor


function __solve(plprob::ProfileLikelihoodProblem, method::FIMProfiler;
  obj0=nothing, parallel_type::Symbol=:none, kwargs...)

  parallel_type == :none || @warn "`parallel_type=$parallel_type` is ignored by `FIMProfiler`; running in serial mode."
  isfinite(plprob.threshold) || throw(ArgumentError("`FIMProfiler` requires finite `plprob.threshold` to define CI level."))

  target = plprob.target
  idxs = get_profile_idxs(target)
  θ= plprob.optpars
  T = float(eltype(θ))
  obj0_t = isnothing(obj0) ? evaluate_obj(plprob, θ) : obj0
  obj_level = obj0_t + plprob.threshold

  F = evaluate_FIM(plprob, θ)
  Fsym = Matrix(Symmetric(F))

  Σ, inversion_used = _invert_matrix(Fsym, get_inversion(method))
  Σ .*= T(get_covariance_factor(method))
  inversion_used != get_inversion(method) &&
    @warn "Requested inversion $(get_inversion(method)) was unstable; used $inversion_used instead."

  z = T(sqrt(plprob.threshold))
  elapsed_time = @elapsed begin
    profile_data = map(idxs) do idx
      lb = get_idx_profile_lb(target, idx)
      ub = get_idx_profile_ub(target, idx)
      θi = θ[idx]
      Σᵢᵢ = T(Σ[idx, idx])

      if !isfinite(Σᵢᵢ)
        x = T[θi]
        obj = T[obj0_t]
        pars = [copy(θ)]
        sol = solution_init(plprob, idx, 0, copy(θ), θi, obj0_t, obj_level)
        return remake(sol; pars=pars, x=x, obj=obj,
          retcodes=(left=:Failure, right=:Failure),
          endpoints=(left=nothing, right=nothing),
          stats=(left=nothing, right=nothing))
      end

      se = sqrt(max(Σᵢᵢ, zero(T)))
      Δ = z * se

      l_raw = θi - Δ
      r_raw = θi + Δ
      l = get_clamp_to_bounds(method) ? clamp(l_raw, lb, ub) : l_raw
      r = get_clamp_to_bounds(method) ? clamp(r_raw, lb, ub) : r_raw

      left_status = (l != l_raw) ? :NonIdentifiable : :Identifiable
      right_status = (r != r_raw) ? :NonIdentifiable : :Identifiable

      x = T[l, θi, r]
      denom = max(Σᵢᵢ, sqrt(eps(T)))
      obj = T[obj0_t + ((xx - θi)^2) / denom for xx in x]
      pars = [begin
          θp = copy(θ)
          θp[idx] = xx
          θp
        end for xx in x]

      sol = solution_init(plprob, idx, 0, copy(θ), θi, obj0_t, obj_level)
      remake(sol; pars=pars, x=x, obj=obj,
        retcodes=(left=left_status, right=right_status),
        endpoints=(left=l, right=r),
        stats=(left=nothing, right=nothing))
    end
  end
  
  return ProfileLikelihoodSolution{typeof(plprob), typeof(profile_data)}(plprob, profile_data, elapsed_time)
end

"""
    evaluate_FIM(plprob::ProfileLikelihoodProblem, θ=plprob.optpars)

Evaluates the Fisher Information Matrix (FIM) at the given parameter values `θ` (default: `plprob.optpars`) by computing the Hessian of the objective function.
"""
function evaluate_FIM(plprob::ProfileLikelihoodProblem,  θ=plprob.optpars)
  return evaluate_hessf(plprob.optprob, θ)
end

function _invert_matrix(F::AbstractMatrix, inversion::Symbol)
  if inversion == :cholesky
    try
      C = cholesky(Symmetric(F), check=true)
      return inv(C), :cholesky
    catch
      return pinv(F), :pinv
    end
  elseif inversion == :pinv
    return pinv(F), :pinv
  #=
  elseif inversion == :svd
    S = svd(F)
    tol = maximum(size(F)) * eps(real(eltype(F))) * maximum(S.S)
    Sinv = Diagonal(map(s -> s > tol ? inv(s) : zero(s), S.S))
    return S.V * Sinv * S.U', :svd
  =#
  end
  throw(ArgumentError("Unsupported inversion mode: $inversion"))
end

