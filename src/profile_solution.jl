
############################################### PROFILE SOLUTION #########################################

struct PLSolution{probType,P}
  prob::probType
  profiles::P
  elapsedTime::Float64
end

function build_profile_solution(plprob::PLProblem, prof_data, elapsed_time)
  PLSolution{typeof(plprob),typeof(prof_data)}(plprob, prof_data, elapsed_time)
end

function Base.show(io::IO, mime::MIME"text/plain", pv::PLSolution) 
  #type_color, no_color = SciMLBase.get_colorizers(io)
  println(io, "Profile Likelihood solution. Use `sol[i]` to index the computed profiles.") 
end

Base.getindex(A::PLSolution, i::Int) = A.profiles[i]
Base.size(A::PLSolution) = size(A.profiles)
Base.length(A::PLSolution) = length(A.profiles)

################################################ PROFILE VALUES ##########################################

# mimics the way saving callback works in SciML DiffEqCallbacks.jl
struct ProfileValues{probType,parsType,xType,objType,retType,epType,statsType}
  dense::Bool
  plprob::probType
  pars::Vector{parsType}
  x::Vector{xType}
  obj::Vector{objType}
  obj_level::Float64
  retcodes::retType
  endpoints::epType
  #=
  retcode_lower::Symbol
  retcode_upper::Symbol
  endpoint_lower::EP1
  endpoint_upper::EP2
  stats_lower::S1
  stats_upper::S2
  =#
  stats::statsType
  #deriv::dobjType
end

#=
set_retcode!(pv::ProfileValues, ::Val{true}, retcode) = pv.lower_retcode[] = retcode
set_retcode!(pv::ProfileValues, ::Val{false}, retcode) = pv.upper_retcode[] = retcode
set_endpoint!(pv::ProfileValues, ::Val{true}, val) = pv.lower_endpoint[] = val
set_endpoint!(pv::ProfileValues, ::Val{false}, val) = pv.upper_endpoint[] = val
=#

get_plprob(pv::ProfileValues) = pv.plprob
get_optprob(pv::ProfileValues) = get_optprob(get_plprob(pv))
get_obj_level(pv::ProfileValues) = pv.obj_level
get_retcodes(pv::ProfileValues) = pv.retcodes
get_endpoints(pv::ProfileValues) = pv.endpoints
get_stats(pv::ProfileValues) = pv.stats

function ProfileValues(::Val{false}, plprob::PLProblem, ::Type{parsType}, ::Type{xType}, ::Type{objType}, obj_level) where {parsType, xType, objType}
  ProfileValues{typeof(plprob), parsType, xType, objType, Nothing, Nothing, Nothing}(false, plprob, Vector{parsType}(), Vector{xType}(), Vector{objType}(), obj_level, nothing, nothing, nothing)
end

function ProfileValues(::Val{false}, plprob::PLProblem, pars::Vector{parsType}, x::Vector{xType}, obj::Vector{objType}, obj_level, retcode, endpoints, stats) where {parsType, xType, objType}
  ProfileValues{typeof(plprob), parsType, xType, objType, typeof(retcode), typeof(endpoints), typeof(stats)}(false, plprob, pars, x, obj, obj_level, retcode, endpoints, stats)
end

isdense(pv::ProfileValues) = pv.dense

function Base.show(io::IO, mime::MIME"text/plain", pv::ProfileValues) 
  #type_color, no_color = SciMLBase.get_colorizers(io)
  println(io, "Profile values. Use `plot` to visualize the profile or `DataFrame` to get tabular representation.") 
  println(io, "Estimated confidence interval (CI): $(get_endpoints(pv))")
  println(io, "CI endpoints retcodes: $(get_retcodes(pv))")
end


function DataFrames.DataFrame(pv::ProfileValues)
  npars = length(pv.pars[1])
  df = DataFrame([getindex.(pv.pars, i) for i in 1:npars], :auto, copycols=false)
  df[!,:objective] = pv.obj
  return df
end
