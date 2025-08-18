"""
   AbstractPLProblem

Abstract type for all profile likelihood problems. Concrete subtypes carry the
information required to define a specific profile likelihood problem.
"""
abstract type AbstractPLProblem{T} end

"""
    AbstractProfileTarget

Abstract type for all profile targets. Concrete subtypes carry the
information required to define what is being profiled. 
"""
abstract type AbstractProfileTarget end

"""
    ParameterTarget{I,V}

Profile target representing profiling of model parameters.

### Fields

- `idxs::AbstractVector{<:Integer}`: Indices of the parameters being profiled.
- `lb::AbstractVector{<:Real}`: Lower bounds for the profile likelihood. 
- `ub::AbstractVector{<:Real}`: Upper bounds for the profile likelihood. 

Profile bounds `lb` and `ub` should be set to finite numerical values.

### Constructors

Create a target with explicit lower and upper bounds for each index.
```julia
ParameterTarget(; idxs::AbstractVector{<:Integer}, lb::AbstractVector{<:Real}, ub::AbstractVector{<:Real})
```
"""
struct ParameterTarget{I,V} <: AbstractProfileTarget
  idxs::I
  lb::V
  ub::V

  function ParameterTarget(idxs::I, lb::AbstractVector{<:Real}, ub::AbstractVector{<:Real}) where {I<:AbstractVector{<:Integer}}
    n = length(idxs)
    (n == length(lb) == length(ub) ) ||
      throw(DimensionMismatch("idxs, lb and ub must have the same length"))
    all(lb .<= ub) || throw(ArgumentError("lower bounds must not exceed upper bounds"))
    T = promote_type(eltype(lb), eltype(ub))
    lb_typed = T.(lb)
    ub_typed = T.(ub)
    new{I,typeof(lb_typed)}(idxs, lb_typed, ub_typed)
  end
end

function Base.show(io::IO, mime::MIME"text/plain", target::ParameterTarget)
  println(io, "ParameterTarget: $(length(target.idxs)) parameters to profile.")
end

############################### SELECTORS ###############################

get_profile_range(prob::ParameterTarget) = tuple.(prob.lb, prob.ub)

############################### CONSTRUCTORS ###############################

# Vector of bounds
ParameterTarget(; idxs::AbstractVector{<:Integer}, lb::AbstractVector{<:Real}, ub::AbstractVector{<:Real}) =
  ParameterTarget(idxs, lb, ub)

#=
# Single parameter
ParameterTarget(idx::Integer; lb::Real, ub::Real) =
  ParameterTarget(; idxs=[idx], lb=[lb], ub=[ub])


# Vector of ranges, one per index
ParameterTarget(idxs::AbstractVector{<:Integer}, profile_range::AbstractVector{<:Tuple{<:Real,<:Real}}) =
  ParameterTarget(idxs, map(first, profile_range), map(last, profile_range))

ParameterTarget(idx::Integer, profile_range::Tuple{<:Real,<:Real}) =
  ParameterTarget([idx], [first(profile_range)], [last(profile_range)])

# Single range applied to all indices
ParameterTarget(idxs::AbstractVector{<:Integer}, range::Tuple{<:Real,<:Real}) =
  ParameterTarget(idxs, fill(first(range), length(idxs)), fill(last(range), length(idxs)))

# Scalar bounds applied to all indices
ParameterTarget(idxs::AbstractVector{<:Integer}, lb::Real, ub::Real) =
  ParameterTarget(idxs, fill(lb, length(idxs)), fill(ub, length(idxs)))
=#


"""
    PLProblem{T,probType,P}
    
Defines a profile likelihood problem.

## Mathematical Specification of a Profile Likelihood Problem:

A profile likelihood problem is defined by 
- an objective function (usually negative log-likelihood function) wrapped within an `optprob::OptimizationProblem`. Consult [Optimization.jl docs](https://docs.sciml.ai/Optimization/stable/API/optimization_problem/) for details.
- a set of optimal values of the parameters `optpars` that minimize the objective function.

### Constructors

```julia
PLProblem(optprob::OptimizationProblem, optpars::AbstractVector{<:Real}, target::AbstractProfileTarget; 
  conf_level::Float64 = 0.95, df::Int = 1, threshold::Real = chi2_quantile(conf_level, df))
```

### Arguments

- `optprob::OptimizationProblem`: The `OptimizationProblem` to be solved.
- `optpars::AbstractVector{<:Real}`: Initial (optimal) values of the parameters.
- `target::AbstractProfileTarget`: The function of parameters to be profiled. Consult the documentation for details.

### Keyword arguments

- `conf_level::Float64`: The confidence level for the profile likelihood. Defaults to `0.95`.
- `df::Int`: The degrees of freedom for the profile likelihood. Defaults to `1`.
- `threshold::Real`: The threshold for the profile likelihood. Can be set to `Inf` if confidence interval endpoint estimation is not required. Defaults to `chi2_quantile(conf_level, df)`.
"""
struct PLProblem{T,probType,P} <: AbstractPLProblem{T}
  optprob::probType
  optpars::P
  target::T
  threshold::Float64
end

sym_profile_type(plprob::PLProblem) = nameof(typeof(get_profile_target(plprob)))
#sym_profile_type(plprob::PLProblem{<:FunctionProfile}) = :function
#isparprofile(::PLProblem{ipp}) where {ipp} = ipp

function Base.show(io::IO, mime::MIME"text/plain", plprob::PLProblem) 
  profile_type = sym_profile_type(plprob)
  println(io, "Profile Likelihood Problem. Profile type: $profile_type")
  println(io, "Parameters' optimal values: ")
  show(io, mime, get_optpars(plprob))
end

############################### SELECTORS ###############################

get_profile_type(prob::PLProblem) = prob.profile_type
get_profile_target(prob::PLProblem) = prob.target
get_optprob(prob::PLProblem) = prob.optprob
get_optpars(prob::PLProblem) = prob.optpars
get_profile_range(prob::PLProblem) = get_profile_range(get_profile_target(prob))
get_threshold(prob::PLProblem) = prob.threshold

############################### CONSTRUCTORS ###############################

function PLProblem(optprob::OptimizationProblem, optpars::AbstractVector{<:Real}, target::AbstractProfileTarget; 
  conf_level::Float64 = 0.95, df::Int = 1, threshold::Real = chi2_quantile(conf_level, df))

  threshold <= 0. && throw(ArgumentError("`threshold` must be positive definite."))
  check_prob_target(optprob, optpars, target)
  # promote_target(optprob, optpars, target)

  PLProblem(optprob, optpars, target, float(threshold))
end

################################ HELPERS ################################

hasthreshold(prob::PLProblem) = isfinite(prob.threshold)

function check_prob_target(optprob::OptimizationProblem, optpars::AbstractVector, target::ParameterTarget) 
  u = optprob.u0
  !(u isa AbstractVector ) && throw(ArgumentError("Expected `optprob.u0` to be a vector (i.e., of `AbstractVector` type)."))
  !(length(u) == length(optpars)) && throw(DimensionMismatch("`optprob.u0` and `optpars` must have the same length."))
  
  idxs = target.idxs
  checkbounds(optpars, idxs)

  lb = target.lb
  ub = target.ub
  
  for idx in idxs
    !(isfinite(lb[idx]) && isfinite(ub[idx])) && throw(ArgumentError("Profile bounds must be finite numbers."))
    !(lb[idx] <= optpars[idx] <= ub[idx]) && throw(ArgumentError("The initial values provided for profiling parameters must lie within the specified `lb[idx] ≤ x[idx] ≤ ub[idx]`"))
  end

  return nothing
end

################################ REMAKE ################################

function SciMLBase.remake(plprob::PLProblem;
  optprob = missing, optpars = missing, target = missing, threshold = missing)
  
  if optprob === missing
    optprob = get_optprob(plprob)
  end
  if optpars === missing
    optpars = get_optpars(plprob)
  end
  if target === missing
    target = get_profile_target(plprob)
  end
  if threshold === missing
    threshold = get_threshold(plprob)
  end

  return PLProblem(optprob, optpars, profile_range; threshold)
end
