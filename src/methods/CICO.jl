Base.@kwdef struct CICOProfiler2{opType, optsType} <: AbstractProfilerMethod
  optimizer::opType
  optimizer_opts::optsType = NamedTuple()
  auglag_opts = NamedTuple()
  resolution::Int = 1
end

mutable struct CICOProfiler2SolverCache{O,StatsType,R} <: AbstractSolverCache
    opt_cache::O
    stats::StatsType
    retcode_original::R
end

function solver_cache_init(plprob::ProfileLikelihoodProblem, target::ParameterTarget, method::CICOProfiler2, idx, dir, profile_range)
    # Implement the initialization logic for the solver cache here
end

function profiler_step!(prof_cache::CICOProfiler2Cache, target::ParameterTarget, solver_cache::CICOProfiler2SolverCache)
    # Implement the step logic for the CICO profiler here
end