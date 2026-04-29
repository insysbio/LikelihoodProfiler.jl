############################## FIMProfiler ##############################

"""
    FIMProfiler

Fisher Information Matrix (FIM)-based asymptotic confidence intervals (Wald approximation).
By default this method reuses Hessian logic from `OptimizationProblem` (user-supplied Hessian or AD backend).
The confidence interval is computed as `θ̂ ± z * sqrt(Σ[idx, idx])`, where 
        - `θ̂` is the `optpars[idx]`, 
        - `z` is the quantile of the chi-squared distribution corresponding to the `conf_level` and `df` parameters of the `ProfileLikelihoodProblem`,
        - `Σ` is the covariance matrix obtained by inverting the FIM.

### Fields

- `inversion::Symbol`: Matrix inversion strategy (`:cholesky`, `:pinv`).
- `clamp_to_bounds::Bool`: Clip estimated interval endpoints to profile bounds.
"""
Base.@kwdef struct FIMProfiler <: AbstractProfilerMethod
  inversion::Symbol = :cholesky
  clamp_to_bounds::Bool = true

  function FIMProfiler(inversion::Symbol, clamp_to_bounds::Bool)
    inversion in (:cholesky, :pinv) ||
      throw(ArgumentError("`inversion` must be one of :cholesky, :pinv (got $inversion)."))

    new(inversion, clamp_to_bounds)
  end
end

get_inversion(fp::FIMProfiler) = fp.inversion
get_clamp_to_bounds(fp::FIMProfiler) = fp.clamp_to_bounds


function __solve(plprob::ProfileLikelihoodProblem, method::FIMProfiler;
  obj0=nothing, parallel_type::Symbol=:none, kwargs...)

  parallel_type == :none || @warn "`parallel_type=$parallel_type` is ignored by `FIMProfiler`; running in serial mode."
  isfinite(plprob.threshold) || throw(ArgumentError("`FIMProfiler` requires finite `plprob.threshold` to define CI level."))

  target = plprob.target
  idxs = get_profile_idxs(target)
  θ̂ = plprob.optpars
  T = float(eltype(θ̂))
  obj0_t = isnothing(obj0) ? evaluate_obj(plprob, θ̂ ): obj0
  obj_level = obj0_t + plprob.threshold

  F = evaluate_FIM(plprob, θ̂)
  Fsym = Matrix(Symmetric(F))

  Σ, inversion_used = _invert_matrix(Fsym, get_inversion(method))
  inversion_used != get_inversion(method) &&
    @warn "Requested inversion $(get_inversion(method)) was unstable; used $inversion_used instead."

  z = T(sqrt(plprob.threshold))
  elapsed_time = @elapsed begin
    profile_data = map(idxs) do idx
      lb = get_idx_profile_lb(target, idx)
      ub = get_idx_profile_ub(target, idx)
      θi = θ̂[idx]
      Σᵢᵢ = T(Σ[idx, idx])

      if !isfinite(Σᵢᵢ)
        x = T[θi]
        obj = T[obj0_t]
        pars = [copy(θ̂)]
        base = solution_init(plprob, idx, 0, copy(θ̂), θi, obj0_t, obj_level)
        return remake(base; pars=pars, x=x, obj=obj,
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
      obj = T[obj0_t + ((xx - θi)^2) / (2 * denom) for xx in x]
      pars = [begin
          θp = copy(θ̂)
          θp[idx] = xx
          θp
        end for xx in x]

      base = solution_init(plprob, idx, 0, copy(θ̂), θi, obj0_t, obj_level)
      remake(base; pars=pars, x=x, obj=obj,
        retcodes=(left=left_status, right=right_status),
        endpoints=(left=l, right=r),
        stats=(left=nothing, right=nothing))
    end
  end

  return ProfileLikelihoodSolution{typeof(plprob), typeof(profile_data)}(plprob, profile_data, elapsed_time)
end

function evaluate_FIM(plprob::ProfileLikelihoodProblem, method::FIMProfiler, θ̂=plprob.optpars)
  return evaluate_hessf(plprob.optprob, θ̂)
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

