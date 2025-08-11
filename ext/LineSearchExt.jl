module LineSearchExt

using LikelihoodProfiler, LineSearch

function compute_step_size(profiler_state::LikelihoodProfiler.ProfilerState, ls_alg::LineSearch.AbstractLineSearchAlgorithm, θ_dir)
  θ_cur = get_curpars(profiler_state)
  opt_prob = get_optprob(profiler_state)
  obj_func = opt_prob.f.f
  fθ_cur = evaluate_optf(opt_prob, θ_cur)

  nl_prob = NonlinearProblem(obj_func, θ_cur)
  ls_cache = init(nl_prob, ls_alg, fθ_cur, θ_cur)

  ls_sol = solve!(ls_cache, θ_cur, θ_dir)

  return ls_sol.step_size, ls_sol.retcode
end

end# LineSearchExt module

