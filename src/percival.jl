# evaluate lower/upper bound of scan_func
function get_endpoint(
  theta_init::Vector{Float64}, # initial parameters' values
  theta_num::Int, # number of parameter to scan

  loss_func::Function, # lambda(theta) - labmbda_min - delta_lambda

  method::Val{:CICO_PERCIVAL},
  direction::Symbol = :upper; 

  theta_bounds::Vector{Tuple{Float64,Float64}} = fill(
      (-Inf, Inf), length(theta_init)
      ),
  scan_bound::Float64 = -9.,

  atol::Float64 = 1e-6,
  rtol::Float64 = 1e-6,
  ctol::Float64 = 1e-6, 

  max_eval::Int = 100000,
  max_time::Float64 = 30.0,
  max_iter::Int = 2000
  #kwargs...
)
  # checking arguments
  if theta_num > length(theta_init)
    throw(DomainError(theta_num, "theta_num exceed theta dimension"))
  end

  islower = direction == :lower
  theta_sign = islower ? 1 : -1
  scan_func(theta) = theta_sign*theta[theta_num]

  # box constrains
  lvar = [tb[1] for tb in theta_bounds]
  uvar = [tb[2] for tb in theta_bounds]

  # nlp problem
  # ADNLPModel(f, x0, lvar, uvar, c, lcon, ucon)
  nlp = ADNLPModel(scan_func, theta_init, lvar, uvar, loss_func, [-Inf], [0.0])


  # callback
  out_of_bound = false
  function cb(nlp, solver, stats)
    if (loss_func(solver.x) < 0.) && (scan_func(solver.x) < scan_bound)
      out_of_bound = true
      stats.status = :user
    end
  end

  # optimizer
  opt = percival(nlp; verbose=0, callback=cb, 
    atol, rtol, ctol, max_eval, max_time, max_iter)

  # optimization results
  optf = theta_sign*opt.objective
  optx = opt.solution
  ret  = opt.status

  if ret == :user && out_of_bound
      pp = ProfilePoint[]
      res = (nothing, pp, :SCAN_BOUND_REACHED)
  elseif ret == :first_order 
      loss = loss_func(optx)
      pp = [ ProfilePoint(optf, loss, optx, ret, nothing) ]
      res = (optf, pp, :BORDER_FOUND_BY_SCAN_TOL)
  else
      pp = ProfilePoint[]
      res = (nothing, pp, ret)
  end

  return res
end
