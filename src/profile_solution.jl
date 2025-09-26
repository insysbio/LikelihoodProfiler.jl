
################################################ PROFILE CURVE ##########################################

const BranchRetcodes = NamedTuple{(:left,:right),
  Tuple{Symbol,Symbol}}

const BranchEndpoints{E} = NamedTuple{(:left,:right),
  Tuple{Union{Nothing, E}, Union{Nothing, E}}}

const BranchStats{S} = NamedTuple{(:left,:right),
  Tuple{Union{Nothing, S}, Union{Nothing, S}}}

# mimics the way saving callback works in SciML DiffEqCallbacks.jl
struct ProfileCurve{P,I,Θ,T,S}
  dense::Bool
  plprob::P
  idx::I
  dir::Int
  pars::Vector{Θ}
  x::Vector{T}
  obj::Vector{T}
  obj_level::T
  retcodes::BranchRetcodes
  endpoints::BranchEndpoints{T}
  stats::BranchStats{S}
end

function SciMLBase.remake(c::ProfileCurve{P,I,Θ,T,S}; 
  dense      = c.dense,
  plprob     = c.plprob,
  idx        = c.idx,
  dir        = c.dir,
  pars       = c.pars,
  x          = c.x,
  obj        = c.obj,
  obj_level  = c.obj_level,
  retcodes   = c.retcodes,
  endpoints  = c.endpoints,
  stats      = c.stats,
) where {P,I,Θ,T,S}
  return ProfileCurve{P,I,Θ,T,S}(dense, plprob, idx, dir, pars, x, obj, obj_level, retcodes, endpoints, stats)
end

solution_init(plprob::ProfileLikelihoodProblem, idx, dir::Int, θ0, x0, obj0, obj_level::T) where {T} = 
  ProfileCurve{typeof(plprob), typeof(idx), typeof(θ0), T, SciMLBase.OptimizationStats}(false, plprob, idx, dir, 
    [θ0], [x0], [obj0], obj_level, (left=:NotStarted, right=:NotStarted), (left=nothing, right=nothing), (left=nothing, right=nothing))

function profiler_init(sol::ProfileCurve, solver_cache::AbstractSolverCache, 
  θ0, x0, obj0, profile_bound, maxiters::Int, verbose::Bool)
  
  stats = solver_cache.stats
  return ProfilerCache{typeof(sol), typeof(θ0), typeof(solver_cache), typeof(x0), typeof(stats)}(
    sol, solver_cache, similar(θ0), copy(θ0), x0, obj0, profile_bound, 1, maxiters, verbose, stats, :Default)
end

isdense(pc::ProfileCurve) = pc.dense
get_plprob(pc::ProfileCurve) = pc.plprob
get_optprob(pc::ProfileCurve) = get_optprob(get_plprob(pc))
get_pars_prev(pc::ProfileCurve) = pc.pars[end-1]
get_x_prev(pc::ProfileCurve) = pc.x[end-1]
get_obj_prev(pc::ProfileCurve) = pc.obj[end-1]

retcodes(pc::ProfileCurve) = pc.retcodes
endpoints(pc::ProfileCurve) = (left = pc.endpoints.left, right = pc.endpoints.right)
stats(pc::ProfileCurve) = (left = pc.stats.left, right = pc.stats.right)
obj_level(pc::ProfileCurve) = pc.obj_level

function Base.show(io::IO, mime::MIME"text/plain", pc::ProfileCurve) 
  #type_color, no_color = SciMLBase.get_colorizers(io)
  branch = if pc.dir == -1 
      " (left-branch)" 
    elseif pc.dir == 1 
      " (right-branch)" 
    else
      ""
    end
  println(io, "Profile curve$(branch). Use `plot` to visualize the profile or `DataFrame` to get tabular representation.") 
  println(io, "Estimated confidence interval (CI): $(endpoints(pc))")
  println(io, "CI endpoints retcodes: $(retcodes(pc))")
end


function DataFrames.DataFrame(pc::ProfileCurve)
  npars = length(pc.pars[1])
  df = DataFrame([getindex.(pc.pars, i) for i in 1:npars], :auto, copycols=false)
  df[!,:objective] = pc.obj
  return df
end

function profiler_save_values!(sol::ProfileCurve, profiler_cache::ProfilerCache)
  @unpack iter, pars_cur, x_cur, obj_cur = profiler_cache

  SciMLBase.copyat_or_push!(sol.pars, iter, pars_cur)
  SciMLBase.copyat_or_push!(sol.x,    iter, x_cur)
  SciMLBase.copyat_or_push!(sol.obj,  iter, obj_cur)
  
  return nothing
end

function reverse_profile!(sol::ProfileCurve)
  reverse!(sol.pars)
  reverse!(sol.x)
  reverse!(sol.obj)
  return nothing
end

function interpolate_endpoint(sol::ProfileCurve)
  # future check if dense
  ol = obj_level(sol)
  interp = LinearInterpolation(sol.x[end-1:end], sol.obj[end-1:end])
  return interp(ol)
end

function profiler_finalize_solution!(profiler_cache::ProfilerCache, sol::ProfileCurve)

  ep = (isidentifiable(profiler_cache) && length(sol.x) ≥ 2) ? float(interpolate_endpoint(sol)) : nothing
  retcode = profiler_cache.retcode
  stats = profiler_cache.stats    

  left_branch = isleft(profiler_cache) 
  left_branch && reverse_profile!(sol)

  new_retcodes = (
    left  = left_branch ? retcode : :NotStarted,
    right = left_branch ? :NotStarted : retcode,
  )
  new_endpoints = (
    left  = left_branch ? ep : nothing,
    right = left_branch ? nothing : ep,
  )
  new_stats = (
    left  = left_branch ? stats : nothing,
    right = left_branch ? nothing : stats,
  )
  profiler_cache.sol = remake(sol; retcodes = new_retcodes, endpoints = new_endpoints, stats = new_stats)
  return nothing
end

function merge_profile_curves(left_profile::ProfileCurve, right_profile::ProfileCurve)

  left_profile.idx == right_profile.idx || throw(ArgumentError("Cannot merge curves for different indices ($(left_profile.idx)) vs ($(right_profile.idx))."))
  (left_profile.dir == -1 && right_profile.dir == 1) ||
        throw(ArgumentError("Expected one left (dir = -1) and one right (dir = +1) branch; got $(left_profile.dir), $(right_profile.dir)."))

  plprob = left_profile.plprob
  obj_level = left_profile.obj_level

  pars = [left_profile.pars; right_profile.pars]
  x = [left_profile.x; right_profile.x]
  obj = [left_profile.obj; right_profile.obj]

  retcodes  = (left  = left_profile.retcodes.left,
                right = right_profile.retcodes.right)
  endpoints = (left  = left_profile.endpoints.left,
                right = right_profile.endpoints.right)
  stats     = (left  = left_profile.stats.left,
                right = right_profile.stats.right)


  return remake(left_profile; dense=false, left_profile.idx, dir=0,  pars=pars, x=x, obj=obj, retcodes=retcodes, endpoints=endpoints, stats=stats)
end


############################################### PROFILE SOLUTION #########################################
"""
    ProfileLikelihoodSolution{probType,P}

Contains the results of a profile likelihood analysis.

### Fields

- `prob::probType`: The profile likelihood problem `ProfileLikelihoodProblem`.
- `profiles::P`: The computed profile curves.
- `elapsed_time::Float64`: The time elapsed during the computation.

### Selectors

A number of selectors are available to extract information from the `sol::ProfileLikelihoodSolution` object. These can be applied to each computed profile `sol[i]`:

- `endpoints(sol[i])`: Returns the confidence interval (CI) endpoints, marking the intersection of the profile with the `threshold`.
- `retcodes(sol[i])`: Returns the retcodes of the CI endpoints estimation.
- `stats(sol[i])`: Returns the statistics of the profile computation.
"""
struct ProfileLikelihoodSolution{probType,P}
  prob::probType
  profiles::P
  elapsed_time::Float64
end

function Base.show(io::IO, mime::MIME"text/plain", pc::ProfileLikelihoodSolution) 
  #type_color, no_color = SciMLBase.get_colorizers(io)
  println(io, "Profile Likelihood solution. Use `sol[i]` to index the computed profiles.") 
end

Base.getindex(A::ProfileLikelihoodSolution, i::Int) = A.profiles[i]
Base.size(A::ProfileLikelihoodSolution) = size(A.profiles)
Base.length(A::ProfileLikelihoodSolution) = length(A.profiles)

function build_profile_solution(plprob::ProfileLikelihoodProblem,
                                branches::AbstractVector,
                                elapsed_time::Real)
  nb = length(branches)
  nb == 0 && return ProfileLikelihoodSolution(plprob, ProfileCurve[], float(elapsed_time))
  iseven(nb) || throw(ArgumentError("Expected an even number of branch curves (left/right per index); got $nb."))

  np = nb ÷ 2
  profile_data = Vector{eltype(branches)}(undef, np)
  @inbounds for i in 1:np
    left_profile  = branches[2i-1]
    right_profile = branches[2i]
    profile_data[i] = merge_profile_curves(left_profile, right_profile)
  end

  ProfileLikelihoodSolution{typeof(plprob),typeof(profile_data)}(plprob, profile_data, float(elapsed_time))
end
