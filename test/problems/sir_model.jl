# ODE function
function ode_func(du,u,p,t)
  S = u[1] 
  I = u[2] 
  R = u[3] 
  b = p[1] 
  g = p[2]
   
  du[1] = -b*S*I
  du[2] = b*S*I - g*I
  du[3] = g*I
end

# we use times, data, parameters and initial values from
# https://github.com/marisae/param-estimation-SIR/blob/master/R/SIR_Example_Main.R#L11-L46
# to build an ODE Problem
times = [0., 7., 14., 21., 28., 35., 42., 49., 56., 63., 70., 77., 84., 91., 98.]
data = [97., 271., 860., 1995., 4419., 6549., 6321., 4763., 2571., 1385., 615., 302., 159., 72., 34.]
tspan = (0., 98.)
p0 = [0.4, 0.25, 1/80000]  

u0_func(p,t) = [1-(data[1]*p[3]), data[1]*p[3], 0.]

ode_prob = ODEProblem(ode_func, u0_func, tspan, p0);

# solver algorithm, tolerances
solver_opts = Dict(
    :alg => AutoTsit5(Rosenbrock23()),
    :reltol => 1e-6,
    :abstol => 1e-8
)

function sir_obj(
  parscur, p;
  ode_prob=ode_prob,
  solver_opts=solver_opts,
  times=times,
  data=data
) 
  # update prob with current parameter values
  probcur = remake(ode_prob,p=parscur)
  
  # solve odes
  sol = solve(
      probcur, 
      solver_opts[:alg]; 
      reltol=solver_opts[:reltol],
      abstol=solver_opts[:abstol],
      saveat=times,
      save_idxs = 2
  )
  
  # current parameter values
  params = probcur.p
  
  # observable
  y = abs.(sol/params[3])
  
  # loss
  return sum(y) - sum(data.*log.(y))
end
