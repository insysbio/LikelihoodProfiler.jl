abstract type AbstractProfilerStep{S} end

const DEFAULT_INIT_STEP = 1e-2

struct FixedStep{S} <: AbstractProfilerStep{S}
  initial_step::S
end 

FixedStep(; initial_step) = FixedStep(prepare_initial_step(initial_step))

struct AdaptiveStep{S} <: AbstractProfilerStep{S}
  reg_order::Int
  initial_step::S
  min_x_step::Float64
  max_x_step::Float64
  min_obj_step::Float64
  max_obj_step::Float64
end 

prepare_initial_step(step::Number) = (p0,i) -> float(step)
prepare_initial_step(step::Function) = step

get_step(s::AbstractProfilerStep{S}, pars, i) where S <: Function = s.initial_step(pars, i)


function propose_next_pars!(profiler::ProfilerState, s::FixedStep)
  paridx = get_idx(profiler)
  curpars = get_curpars(profiler)
  dir = get_dir(profiler)
  step = get_step(s, curpars, paridx)
  bound = get_profile_bound(profiler)

  profiler.pars_cache .= curpars
  profiler.pars_cache[paridx] = crop_to_bound(curpars[paridx] + dir*step, dir, bound)
  return profiler.pars_cache
end

function propose_next_pars!(profiler::ProfilerState, s::AdaptiveStep)
  #TODO: Implement adaptive step 
end

propose_next_pars!(profiler::ProfilerState, s::AdaptiveStep, reg_order::Val{0}) = 
  propose_next_pars!(profiler, FixedStep(s.initial_step))


##################### HELPER FUNCTIONS #####################

function crop_to_bound(xguess::Number, dir::Int, bound::Number)
  if dir*xguess > dir*bound
    return bound
  else
    return xguess
  end
end