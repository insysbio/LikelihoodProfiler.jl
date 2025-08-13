profile(plprob::PLProblem{ParameterProfile}, method::AbstractProfilerMethod; 
  idxs::AbstractVector{<:Int} = eachindex(get_optpars(plprob)),
  parallel_type::Symbol=:none, kwargs...) = solve(plprob, method; idxs=idxs, parallel_type=parallel_type, kwargs...)

Base.depwarn("`profile` is deprecated. Use `solve` from CommonSolve.jl instead.", :profile)

