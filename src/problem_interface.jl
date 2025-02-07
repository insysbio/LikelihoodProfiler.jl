
abstract type AbstractPLProblem{T,PF} end

"""

Defines a profile likelihood problem.

## Mathematical Specification of a Profile Likelihood Problem:

A profile likelihood problem is defined by 
- an objective function (usually negative log-likelihood function) wrapped within an `optprob::OptimizationProblem` (see Optimization.jl docs https://docs.sciml.ai/Optimization/stable/API/optimization_problem/).
- a set of parameters `optpars` that minimize the objective.

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
- `threshold`: The threshold for the profile likelihood. Defaults to `chi2_quantile(conf_level, df)`.
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
#get_direction(prob::PLProblem) = prob.direction

# moved to profile() level
#get_conf_level(prob::PLProblem) = prob.conf_level
#get_df(prob::PLProblem) = prob.df
#get_optobj(prob::PLProblem) = prob.optobj
#get_obj_level(prob::PLProblem) = prob.obj_level

############################### CONSTRUCTORS ###############################

function PLProblem(optprob::OptimizationProblem, optpars::AbstractVector{<:Real},
  profile_range = tuple.(optprob.lb, optprob.ub); 
  conf_level::Float64 = 0.95, df::Int = 1, threshold::Real = chi2_quantile(conf_level, df))

  threshold <= 0. && throw(ArgumentError("`threshold` must be positive definite."))
  numpars = length(optpars)
  validate_optpars(numpars, optprob.u0)
  #TODO check if opt_range is not inside profile_range
  Tprofile_range = promote_profile_range(optpars, profile_range)

  build_plproblem(ParameterProfile(), optprob, optpars, nothing, Tprofile_range, float(threshold))
end

#=
# template for FunctionProfile
function PLProblem(optprob::OptimizationProblem, optpars::AbstractVector{<:Real},
  parfunc::OptimizationFunction, profile_range::Tuple{Real,Real};
  conf_level::Float64 = 0.95, df::Int = 1, threshold::Real = chi2_quantile(conf_level, df))

  threshold <= 0. && throw(ArgumentError("`threshold` must be positive definite."))
  numpars = length(optpars)
  validate_optpars(numpars, optprob.u0)
  Tprofile_range = promote_profile_range(compute_optf(parfunc, optpars), profile_range)

  build_plproblem(FunctionProfile(), optprob, optpars, parfunc, Tprofile_range, float(threshold))
end
=#

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

function promote_profile_range(x::AbstractVector, profile_range::AbstractVector)
  length(profile_range) != length(x) && 
    throw(DimensionMismatch("`profile_range` must be either Tuple or AbstractVector of the same size as `optpars`."))
  promote_profile_range.(x, profile_range)
end

promote_profile_range(x::AbstractVector, profile_range::Tuple) = promote_profile_range.(x, Ref(profile_range))

function promote_profile_range(x::Number, profile_range::Tuple)
  lb, ub = profile_range
  validate_profile_bound(lb)
  validate_profile_bound(ub)
  !(lb <= x <= ub) && throw(ArgumentError("Inconsistent `profile_range`: profiling value `x` must be lb <= x <= ub ."))
  return (float(lb),float(ub))
end

function validate_profile_bound(bound) 
  (isnothing(bound) || isinf(bound)) &&
    throw(ArgumentError("`profile_range` must contain finite values."))
end

function validate_optpars(npars::Int, u)
  npars < 2 && throw(DimensionMismatch("`optpars` length must be ≥ 2."))
  !(u isa AbstractVector && npars == length(u)) && 
    throw(ArgumentError("OptimizationProblem initial values must be of the same type and size as `optpars`."))
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

#=
function SciMLBase.remake(plprob::PLProblem{<:FunctionProfile};
  optprob = missing, optpars = missing, parfunc = missing, profile_range = missing, threshold = missing)

  if optprob === missing
    optprob = get_optprob(plprob)
  end
  if optpars === missing
    optpars = get_optpars(plprob)
  end
  if parfunc === missing
    parfunc = get_parfunc(plprob)
  end
  if profile_range === missing
    profile_range = get_profile_range(plprob)
  end
  if threshold === missing
    threshold = get_threshold(plprob)
  end

  return PLProblem(optprob, optpars, parfunc, profile_range; threshold)
end
=#