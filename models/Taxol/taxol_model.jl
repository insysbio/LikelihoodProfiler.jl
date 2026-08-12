
# https://github.com/marisae/cancer-chemo-identifiability/blob/master/Profile%20Likelihood/testa0_de.m
using Distributions

const Vt      = 10.515*100
const V0      = 1.3907*Vt
const lam     = 9.5722
const aRP     = 20
const arstexp = 3
const adthexp = 4
const theta   = 10

function taxol_ode(du, u, p, t)
  a0, ka, r0, d0, kd = p.params
  drug = p.drug

  Ncel = u[:P] + u[:Ap] + u[:Ar]
  Lfac = ((Vt-Ncel)^theta)/((V0^theta) + ((Vt-Ncel)^theta))

  arst = a0*(drug^arstexp)/(ka^arstexp + (drug^arstexp))
  adth = d0*(drug^adthexp)/(kd^adthexp + (drug^adthexp))
  arcv = r0

  du[:P] = -lam*u[:P] + aRP*u[:Ap]*Lfac - arst*u[:P] + arcv*u[:Ar]
  du[:Ap] = 2*lam*u[:P] - aRP*u[:Ap]*Lfac
  du[:Ar] = arst*u[:P] - adth*u[:Ar] - arcv*u[:Ar]
  
  return nothing
end

u0 = ComponentArray(
  P = 7.2700, 
  Ap = 2.5490, 
  Ar = 0.0
)

function taxol_params(x, d)
  return ComponentArray(
    params = ComponentArray(
      a0 = x[1],
      ka = x[2],
      r0 = x[3],
      d0 = x[4],
      kd = x[5]),
    drug = d
  )
end

p0 = taxol_params([8.3170, 8.0959, 0.0582, 1.3307, 119.1363], 5.0)

tspan = (0.,15.)
ode_prob = SciMLBase.ODEProblem(taxol_ode, u0, tspan, p0)

# initial values and parameters
# https://github.com/marisae/cancer-chemo-identifiability/blob/master/Profile%20Likelihood/testa0_soln.m#L3-L6
# https://github.com/marisae/cancer-chemo-identifiability/blob/master/Profile%20Likelihood/testa0_fit.m#L4

#P0 = 7.2700
#R0 = 2.5490

u0 = [7.2700, 2.5490, 0.]
p0 = log10.([8.3170, 8.0959, 0.0582, 1.3307, 119.1363])

tspan = (0.,15.)

prob = SciMLBase.ODEProblem((du,u,p,t)->ode_func(du,u,p,t,5.0), u0, tspan, p0)
 
times = [0., 3., 6., 9., 12., 15.]   # days
dose = [5., 10., 40., 100.];    # dose in ng/ml

# Control data
Cell = [0.009, 0.050, 0.120, 0.189, 0.230, 0.260]*1091.0   # thousands of cells
Cerr = [0.006, 0.012, 0.010, 0.011, 0.011, 0.011]*1091.0   # thousands of cells

# 0.005 ug/ml Taxol
Cell005 = [0.009, 0.047, 0.089, 0.149, 0.198, 0.219]*1091.0   # thousands of cells
Cerr005 = [0.006, 0.013, 0.010, 0.011, 0.013, 0.010]*1091.0   # thousands of cells

# 0.010 ug/ml Taxol
Cell010 = [0.009, 0.043, 0.077, 0.093, 0.109, 0.128]*1091.0   # thousands of cells
Cerr010 = [0.006, 0.012, 0.013, 0.012, 0.014, 0.012]*1091.0   # thousands of cells

# 0.040 ug/ml Taxol
Cell040 = [0.009, 0.025, 0.047, 0.054, 0.076, 0.085]*1091.0   # thousands of cells
Cerr040 = [0.005, 0.010, 0.010, 0.011, 0.010, 0.010]*1091.0   # thousands of cells

# 0.100 ug/ml Taxol
Cell100 = [0.009, 0.025, 0.026, 0.028, 0.029, 0.031]*1091.0   # thousands of cells
Cerr100 = [0.006, 0.010, 0.009, 0.008, 0.011, 0.011]*1091.0   # thousands of cells

C005 = mean(Cell005)
C010 = mean(Cell010)
C040 = mean(Cell040)
C100 = mean(Cell100)

data = [Cell005/C005, Cell010/C010, Cell040/C040, Cell100/C100]
datamean = [C005, C010, C040, C100]

solver_opts = (
    alg = Tsit5(),
    reltol = 1e-6,
    abstol = 1e-8,
)

function taxol_obj(x, hp)

  loss = 0.
  for (i,d) in enumerate(dose)
    prob = remake(ode_prob; p = taxol_params(exp10.(x), d))
    sol = solve(
      prob, 
      solver_opts.alg, 
      reltol=solver_opts.reltol,
      abstol=solver_opts.abstol,
      saveat=times)
    if !SciMLBase.successful_retcode(sol)
     return Inf
    end
      
    for time_idx in eachindex(data[i])
      sim = (sol[1, time_idx] + sol[2, time_idx] + sol[3, time_idx]) / datamean[i]
      loss += abs2(sim - data[i][time_idx])
    end
  end
  return loss
end

taxol_relative_errors = vcat(
    Cerr005 / C005,
    Cerr010 / C010,
    Cerr040 / C040,
    Cerr100 / C100,
)
sigmasq = mean(taxol_relative_errors)^2