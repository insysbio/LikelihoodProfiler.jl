abstract type AbstractProfilerStep{S} end

const DEFAULT_INIT_STEP = 1e-2

struct FixedStep{S} <: AbstractProfilerStep{S}
  initial_step::S
end 

FixedStep(; initial_step) = FixedStep(prepare_initial_step(initial_step))

struct AdaptiveStep{S} <: AbstractProfilerStep{S}
  reg_order::Int
  initial_step::S
  max_step::Float64
  min_step::Float64
end 

prepare_initial_step(step::Number) = (p0,i) -> float(step)
prepare_initial_step(step::Function) = step

#get_step(s::AbstractProfilerStep{S}, pars, i) where S <: Number = s.step
get_step(s::AbstractProfilerStep{S}, pars, i) where S <: Function = s.initial_step(pars, i)

function compute_next_pars!(profiler::ProfilerState, s::FixedStep)
  
  paridx = get_idx(profiler)
  curpars = get_curpars(profiler)
  dir = get_dir(profiler)
  step = get_step(s, curpars, paridx)
  bound = get_profile_bound(profiler)

  profiler.parscache .= curpars
  profiler.parscache[paridx] = cut_xguess(curpars[paridx] + dir*step, dir, bound)
  return profiler.parscache
end

function compute_next_pars!(profiler::ProfilerState, s::AdaptiveStep)
  #TODO: Implement adaptive step 
end

##################### HELPER FUNCTIONS #####################

function cut_xguess(xguess::Number, dir::Int, bound::Number)
  if dir*xguess > dir*bound
    return bound
  else
    return xguess
  end
end