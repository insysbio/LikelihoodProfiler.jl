
abstract type AbstractPLProblem{T,PF} end

"""
    PLProblem{T,probType,P,PF,PR}
    
Defines a profile likelihood problem.

## Mathematical Specification of a Profile Likelihood Problem:

A profile likelihood problem is defined by 
- an objective function (usually negative log-likelihood function) wrapped within an `optprob::OptimizationProblem`. Consult [Optimization.jl docs](https://docs.sciml.ai/Optimization/stable/API/optimization_problem/) for details.
- a set of optimal values of the parameters `optpars` that minimize the objective function.

### Constructors

```julia
PLProblem(optprob, optpars, profile_range = tuple.(optprob.lb, optprob.ub); 
  conf_level = 0.95, df = 1, threshold = chi2_quantile(conf_level, df))
```

### Arguments

- `optprob`: The `OptimizationProblem` to be solved.
- `optpars`: Initial (optimal) values of the parameters.
- `profile_range`: The range over which the profile likelihood is computed. Defaults to `tuple.(lb,ub)` of the `OptimizationProblem`.

### Keyword arguments

- `conf_level`: The confidence level for the profile likelihood. Defaults to `0.95`.
- `df`: The degrees of freedom for the profile likelihood. Defaults to `1`.
- `threshold`: The threshold for the profile likelihood. Can be set to `Inf` if confidence interval endpoint estimation is not required. Defaults to `chi2_quantile(conf_level, df)`.
"""
struct PLProblem{T,probType,P,PF,PR} <: AbstractPLProblem{T,PF}
  profile_type::T
  optprob::probType
  optpars::P
  parfunc::PF
  profile_range::PR
#  direction::D
  threshold::Float64
end

sym_profile_type(plprob::PLProblem{<:ParameterProfile}) = :parameter
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
get_optprob(prob::PLProblem) = prob.optprob
get_optpars(prob::PLProblem) = prob.optpars
get_parfunc(prob::PLProblem) = prob.parfunc
get_profile_range(prob::PLProblem) = prob.profile_range
get_threshold(prob::PLProblem) = prob.threshold

############################### CONSTRUCTORS ###############################

function PLProblem(optprob::OptimizationProblem, optpars::AbstractVector{<:Real},
  profile_range::Union{AbstractVector, Tuple} = tuple.(optprob.lb, optprob.ub); 
  conf_level::Float64 = 0.95, df::Int = 1, threshold::Real = chi2_quantile(conf_level, df))

  threshold <= 0. && throw(ArgumentError("`threshold` must be positive definite."))

  validate_dims(optprob.u0, optpars, profile_range)
  #TODO check if opt_range is not inside profile_range
  #Tprofile_range = promote_profile_range(optpars, profile_range)

  build_plproblem(ParameterProfile(), optprob, optpars, nothing, profile_range, float(threshold))
end

function build_plproblem(
  profile_type::T,
  optprob::OptimizationProblem,
  optpars::AbstractVector{<:Real},
  parfunc::Union{Nothing,OptimizationFunction},
  profile_range,
  threshold::Float64
) where T
  return PLProblem{T,typeof(optprob), typeof(optpars), typeof(parfunc), typeof(profile_range)}(
    profile_type, optprob, optpars, parfunc, profile_range, threshold)
end

################################ HELPERS ################################

hasthreshold(prob::PLProblem) = isfinite(prob.threshold)

function validate_dims(u, optpars::AbstractVector, profile_range) 
  !(u isa AbstractVector ) && throw(ArgumentError("Expected `optprob.u0` to be a vector (i.e., of `AbstractVector` type)."))
  !(length(u) == length(optpars)) && throw(DimensionMismatch("`optprob.u0` and `optpars` must have the same length."))
  if profile_range isa AbstractVector
    !(length(profile_range) == length(optpars)) && 
      throw(DimensionMismatch("`profile_range` must be either Tuple or AbstractVector of the same size as `optpars`."))
  end
  return nothing
end

################################ REMAKE ################################

function SciMLBase.remake(plprob::PLProblem{<:ParameterProfile};
  optprob = missing, optpars = missing, profile_range = missing, threshold = missing)
  
  if optprob === missing
    optprob = get_optprob(plprob)
  end
  if optpars === missing
    optpars = get_optpars(plprob)
  end
  if profile_range === missing
    profile_range = get_profile_range(plprob)
  end
  if threshold === missing
    threshold = get_threshold(plprob)
  end

  return PLProblem(optprob, optpars, profile_range; threshold)
end
