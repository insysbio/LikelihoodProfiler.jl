
"""
    ParameterTarget{I,B}

Profile target representing profiling of model parameters.

### Fields

- `idxs::AbstractVector{<:Integer}`: Indices of the parameters being profiled.
- `profile_lower::AbstractVector{<:Real}`: Lower bounds for the profile likelihood. 
- `profile_upper::AbstractVector{<:Real}`: Upper bounds for the profile likelihood. 

Profile bounds `profile_lower` and `profile_upper` should be vectors of finite numerical values. 

### Constructors

Create a target with explicit lower and upper bounds for each index.
```julia
ParameterTarget(; idxs::AbstractVector{<:Integer}, profile_lower::AbstractVector{<:Real}, profile_upper::AbstractVector{<:Real})
```
"""
struct ParameterTarget{I<:AbstractVector{<:Integer}, B<:AbstractVector{<:Real}} <: AbstractProfileTarget
  idxs::I
  profile_lower::B
  profile_upper::B

  function ParameterTarget(idxs::I, 
                           profile_lower::AbstractVector{<:Real}, 
                           profile_upper::AbstractVector{<:Real}) where {I<:AbstractVector{<:Integer}}
    n = length(idxs)
    n > 0           || throw(ArgumentError("`idxs` must be non-empty."))
    allunique(idxs) || throw(ArgumentError("`idxs` must be unique."))
    all(>=(1), idxs) || throw(ArgumentError("`idxs` must be positive integers."))
    
    lb_typed, ub_typed = _promote_profile_bounds(profile_lower, profile_upper, n)
    new{I,typeof(lb_typed)}(idxs, lb_typed, ub_typed)
  end
end

function Base.show(io::IO, mime::MIME"text/plain", target::ParameterTarget)
  println(io, "ParameterTarget: $(length(get_profile_idxs(target))) parameter(s) to profile.")
end

"""
    FunctionTarget{F,B}

Profile target representing profiling of functions of model parameters.

### Fields

- `fs::AbstractVector{<:OptimizationFunction}`: Functions of the parameters being profiled.
- `profile_lower::AbstractVector{<:Real}`: Lower bounds for the profile likelihood. 
- `profile_upper::AbstractVector{<:Real}`: Upper bounds for the profile likelihood. 

Profile bounds `profile_lower` and `profile_upper` should be vectors of finite numerical values. 

### Constructors

Create a target with explicit lower and upper bounds for each function of parameters.
```julia
FunctionTarget(; fs::AbstractVector{<:OptimizationFunction}, profile_lower::AbstractVector{<:Real}, profile_upper::AbstractVector{<:Real})
```
"""
struct FunctionTarget{F<:AbstractVector{<:OptimizationFunction}, B<:AbstractVector{<:Real}} <: AbstractProfileTarget
  fs::F
  profile_lower::B
  profile_upper::B

  function FunctionTarget(fs::F, 
                           profile_lower::AbstractVector{<:Real}, 
                           profile_upper::AbstractVector{<:Real}) where {F<:AbstractVector{<:OptimizationFunction}}
    n = length(fs)
    n > 0 || throw(ArgumentError("`fs` must be non-empty."))
    lb_typed, ub_typed = _promote_profile_bounds(profile_lower, profile_upper, n)
    new{F,typeof(lb_typed)}(fs, lb_typed, ub_typed)
  end
end

function Base.show(io::IO, mime::MIME"text/plain", target::FunctionTarget)
  println(io, "FunctionTarget: $(length(target)) function(s) to profile.")
end
############################### SELECTORS ###############################

get_profile_idxs(t::ParameterTarget) = t.idxs
get_profile_idxs(t::FunctionTarget) = 1:length(t.fs)
Base.length(t::AbstractProfileTarget) = length(get_profile_idxs(t))
get_profile_lb(t::AbstractProfileTarget) = t.profile_lower
get_profile_ub(t::AbstractProfileTarget) = t.profile_upper
get_profile_lb(t::AbstractProfileTarget, idx) = t.profile_lower[idx]
get_profile_ub(t::AbstractProfileTarget, idx) = t.profile_upper[idx]

get_profile_fs(ft::FunctionTarget) = ft.fs

############################### CONSTRUCTORS ###############################

ParameterTarget(; idxs::AbstractVector{<:Integer}, profile_lower::AbstractVector{<:Real}, profile_upper::AbstractVector{<:Real}) =
  ParameterTarget(idxs, profile_lower, profile_upper)

FunctionTarget(; fs::AbstractVector{<:OptimizationFunction}, profile_lower::AbstractVector{<:Real}, profile_upper::AbstractVector{<:Real}) =
  FunctionTarget(fs, profile_lower, profile_upper)

################################ HELPERS ################################

function _promote_profile_bounds(lb::AbstractVector{<:Real}, ub::AbstractVector{<:Real}, n::Int)
  (n == length(lb) == length(ub) ) ||
            throw(DimensionMismatch("`profile_lower` and `profile_upper` must have the same length equal to the number of quantities being profiled."))
  all(isfinite.(lb)) && all(isfinite.(ub)) ||
            throw(ArgumentError("profile bounds must be finite numbers"))
  T = promote_type(float(eltype(lb)), float(eltype(ub)))
  lb_typed = T.(lb); ub_typed = T.(ub)
  all(lb_typed .<= ub_typed) || throw(ArgumentError("`profile_lower` must not exceed `profile_upper`"))
  return lb_typed, ub_typed
end

"""
    ProfileLikelihoodProblem{T,probType,P}
    
Defines a profile likelihood problem.

## Mathematical Specification of a Profile Likelihood Problem:

A problem is specified by:
  • `optprob::OptimizationProblem` — wraps your objective (e.g. negative log-likelihood)
  • `optpars::AbstractVector{<:Real}` — parameter values to start profiling from (typically the optimum)
  • `target::AbstractProfileTarget` — what to profile (parameters or functions)

### Constructors

1. Explicit target interface (advanced)
```julia
ProfileLikelihoodProblem(optprob::OptimizationProblem, optpars::AbstractVector{<:Real}, target::AbstractProfileTarget; 
  conf_level::Float64 = 0.95, df::Int = 1, threshold::Union{Nothing,Real} = nothing)
```
  • `target`: `ParameterTarget` (see [`ParameterTarget`](@ref)) or `FunctionTarget` (see [`FunctionTarget`](@ref)) defining what to profile and the profile bounds.
  • `conf_level`: Confidence level for the profile likelihood. Defaults to `0.95`.
  • `df`: Degrees of freedom for the profile likelihood. Defaults to `1`.
  • `threshold`: Profile likelihood threshold. If not provided, computed from `conf_level` and `df`. Can be set to `Inf` if confidence interval endpoint estimation is not required.

2. Parameter profiling sugar
```julia
ProfileLikelihoodProblem(optprob::OptimizationProblem, optpars::AbstractVector{<:Real};
  idxs = nothing, profile_lower = nothing, profile_upper = nothing, kwargs...)
```
  • `idxs`: Indices of parameters to profile; Integer or vector of integers; if nothing, profile all parameters.
  • `profile_lower`, `profile_upper`: Bounds for profiling. Accept scalars or vectors of finite numbers; if `nothing`, taken from `optprob`. 
    If scalar bounds are provided, they will be expanded to match the number of parameters being profiled.
  • `kwargs...`: passed to the explicit target constructor.

3. Function profiling sugar
```julia
ProfileLikelihoodProblem(optprob::OptimizationProblem, optpars::AbstractVector{<:Real};
  fs = nothing, profile_lower = nothing, profile_upper = nothing, kwargs...)
```
  • `fs`: `OptimizationFunction` or vector of `OptimizationFunction` - functions of parameters to be profiled.
  • `profile_lower`, `profile_upper`: Bounds for profiling. Accept scalars or vectors of finite numbers. 
    If scalar bounds are provided, they will be expanded to match the number of functions being profiled.
  • `kwargs...`: passed to the explicit target constructor.
"""
struct ProfileLikelihoodProblem{T,probType,P} <: AbstractProfileLikelihoodProblem
  optprob::probType
  optpars::P
  target::T
  threshold::Float64
end

function Base.show(io::IO, mime::MIME"text/plain", plprob::ProfileLikelihoodProblem) 
  println(io, "Profile Likelihood Problem. Profile threshold: $(plprob.threshold)")
  show(io, mime, plprob.target)
  println(io, "Parameters' optimal values: ")
  show(io, mime, plprob.optpars)
end

############################### CONSTRUCTORS ###############################

function ProfileLikelihoodProblem(optprob::OptimizationProblem, optpars::AbstractVector{<:Real}, target::AbstractProfileTarget; 
  conf_level::Float64 = 0.95, df::Int = 1, threshold::Union{Nothing, Real} = nothing)

  0 < conf_level < 1 || throw(ArgumentError("`conf_level` must be in (0,1)."))
  df > 0 || throw(ArgumentError("`df` must be a positive integer."))
  if !isnothing(threshold) && !(threshold > 0)
    throw(ArgumentError("`threshold` must be strictly positive. If confidence interval endpoint estimation is not required, set `threshold=Inf`."))
  end

  u = optprob.u0
  !(u isa AbstractVector ) && throw(ArgumentError("Expected `optprob.u0` to be a vector (i.e., of `AbstractVector` type)."))
  !(length(u) == length(optpars)) && throw(DimensionMismatch("`optprob.u0` and `optpars` must have the same length."))

  _check_prob_target(optprob, optpars, target)
  τ = isnothing(threshold) ? chi2_quantile(conf_level, df) : float(threshold)
  ProfileLikelihoodProblem(optprob, optpars, target, τ)
end

function ProfileLikelihoodProblem(optprob::OptimizationProblem, optpars::AbstractVector{<:Real}; 
  idxs=nothing, profile_lower=nothing, profile_upper=nothing, kwargs...)
  
  n = length(optpars)
  lb = isnothing(profile_lower) ? optprob.lb : profile_lower
  ub = isnothing(profile_upper) ? optprob.ub : profile_upper
  (isnothing(lb) || isnothing(ub)) &&
     throw(ArgumentError("Bounds not found in `OptimizationProblem`; pass `profile_lower`/`profile_upper`."))

  I = _materialize_idxs(idxs, n)
  lbv = _materialize_profile_bound(lb, length(I))
  ubv = _materialize_profile_bound(ub, length(I))

  lbI = (length(lbv) == n) ? lbv[I] : lbv
  ubI = (length(ubv) == n) ? ubv[I] : ubv
  lbv_typed, ubv_typed = _ensure_real_bounds(lbI, ubI)
  
  tgt = ParameterTarget(; idxs=I, profile_lower=lbv_typed, profile_upper=ubv_typed)
  return ProfileLikelihoodProblem(optprob, optpars, tgt; kwargs...)
end

function ProfileLikelihoodProblem(optprob::OptimizationProblem, optpars::AbstractVector{<:Real}, fs::Union{F, AbstractVector{F}}; 
  profile_lower=nothing, profile_upper=nothing, kwargs...) where {F<:OptimizationFunction}

  (isnothing(profile_lower) || isnothing(profile_upper)) &&
    throw(ArgumentError("Function targets require `profile_lower` and `profile_upper` (scalar or vector)."))
  Fs = fs isa AbstractVector ? collect(fs) : [fs]

  lbv = _materialize_profile_bound(profile_lower, length(Fs))
  ubv = _materialize_profile_bound(profile_upper, length(Fs))
  lbv_typed, ubv_typed = _ensure_real_bounds(lbv, ubv)

  tgt = FunctionTarget(; fs=Fs, profile_lower=lbv_typed, profile_upper=ubv_typed)
  return ProfileLikelihoodProblem(optprob, optpars, tgt; kwargs...)
end

################################ HELPERS ################################

hasthreshold(prob::ProfileLikelihoodProblem) = isfinite(prob.threshold)

function _check_prob_target(optprob::OptimizationProblem, optpars::AbstractVector, target::ParameterTarget) 

  idxs = get_profile_idxs(target)
  checkbounds(Bool, optpars, idxs) || throw(ArgumentError("`idxs` must be within 1:$(length(optpars))"))

  lb = get_profile_lb(target); ub = get_profile_ub(target)
  for (i, idx) in enumerate(idxs) 
    (lb[i] <= optpars[idx] <= ub[i]) || throw(ArgumentError("Initial value optpars[$idx] = $(optpars[idx]) must satisfy profile_lower[$i] ≤ x ≤ profile_upper[$i]"))
  end
  return nothing
end

function _check_prob_target(optprob::OptimizationProblem, optpars::AbstractVector, target::FunctionTarget) 
  return nothing
end

function _materialize_idxs(idxs, n::Int)
  if isnothing(idxs)
    return collect(Base.OneTo(n))
  elseif idxs isa Integer && 1 <= idxs <= n
    return [Int(idxs)]
  elseif idxs isa AbstractVector{<:Integer} && all(1 .<= idxs .<= n)
    return Int.(collect(idxs))
  else
    throw(ArgumentError("`idxs` must be a positive Integer or a vector of positive Integers."))
  end
end

function _materialize_profile_bound(b, I::Int)
  if b isa Real
    return fill(float(b), I)
  elseif isa(b, AbstractVector)
    return collect(b)
  else
    throw(ArgumentError("`profile_lower` and `profile_upper` must be scalars or vectors of numbers."))
  end
end

function _ensure_real_bounds(lb::AbstractVector, ub::AbstractVector)
  all(x -> x isa Real, lb) && all(x -> x isa Real, ub) ||
    throw(ArgumentError("`profile_lower` and `profile_upper` must be scalars or vectors of Real numbers for selected indices."))
  return float.(lb), float.(ub)
end

################################ REMAKE ################################

function SciMLBase.remake(plprob::ProfileLikelihoodProblem;
  optprob = missing, optpars = missing, target = missing, threshold = missing)
  
  if optprob === missing
    optprob = plprob.optprob
  end
  if optpars === missing
    optpars = plprob.optpars
  end
  if target === missing
    target = plprob.target
  end
  if threshold === missing
    threshold = plprob.threshold
  end

  return ProfileLikelihoodProblem(optprob, optpars, target; threshold)
end
