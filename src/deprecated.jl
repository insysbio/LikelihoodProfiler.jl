function profile(plprob::ProfileLikelihoodProblem{ParameterProfile}, method::AbstractProfilerMethod; 
  idxs::AbstractVector{<:Int} = eachindex(get_optpars(plprob)),
  parallel_type::Symbol=:none, kwargs...) 
  Base.depwarn("`profile` is deprecated. Please use `solve` instead.", :profile, force = true)

  return solve(plprob, method; idxs=idxs, parallel_type=parallel_type, kwargs...)
end


function PLProblem(optprob::OptimizationProblem, optpars::AbstractVector{<:Real},
  profile_range::Union{AbstractVector, Tuple} = tuple.(optprob.lb, optprob.ub); 
  conf_level::Float64 = 0.95, df::Int = 1, threshold::Real = chi2_quantile(conf_level, df))
  Base.depwarn("`PLProblem` is deprecated. Please use `ProfileLikelihoodProblem` instead.", :PLProblem, force = true)

  return ProfileLikelihoodProblem(optprob, optpars, profile_range; 
    conf_level=conf_level, df=df, threshold=threshold)
end
