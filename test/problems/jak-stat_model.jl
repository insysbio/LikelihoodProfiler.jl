
# constants
const cyt = 1.4
const nuc = 0.45
const ratio = 0.693
const specC17 = 0.107


# ode system
function stat5_ode(du, u, p, time)
    # 8 states:
    (STAT5A, pApA, STAT5B, pApB, pBpB, nucpApA, nucpApB, nucpBpB) = u
    # 6 parameters
    (Epo_degradation_BaF3, k_exp_hetero, k_exp_homo, k_imp_hetero, k_imp_homo, k_phos) =  p
    
    BaF3_Epo = 1.25e-7*exp(-1*Epo_degradation_BaF3*time)

    v1 = BaF3_Epo*(STAT5A^2)*k_phos
    v2 = BaF3_Epo*STAT5A*STAT5B*k_phos
    v3 = BaF3_Epo*(STAT5B^2)*k_phos
    v4 = k_imp_homo*pApA
    v5 = k_imp_hetero*pApB
    v6 = k_imp_homo*pBpB
    v7 = k_exp_homo*nucpApA
    v8 = k_exp_hetero*nucpApB
    v9 = k_exp_homo*nucpBpB

    du[1] = -2*v1 - v2 + 2*v7*(nuc/cyt) + v8*(nuc/cyt)
    du[2] = v1 - v4
    du[3] = -v2 -2*v3 + v8*(nuc/cyt) + 2*v9*(nuc/cyt)
    du[4] = v2 - v5
    du[5] = v3 - v6
    du[6] = v4*(cyt/nuc) - v7
    du[7] = v5*(cyt/nuc) - v8
    du[8] = v6*(cyt/nuc) - v9
end;
     
data = DataFrame(
  time = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 160.0, 200.0, 240.0],
  pSTAT5A_rel = [7.901072999, 66.36349397, 81.17132392, 94.73030806, 95.11648305, 91.44171655, 91.25709923, 93.67229784, 88.75423282, 85.26970322, 81.13239534, 76.13592848, 65.24805913, 42.59965871, 25.15779754, 15.4301824],
  pSTAT5B_rel = [4.596533343, 29.63454599, 46.04380647, 81.97473362, 80.5716093, 79.03571964, 75.67238037, 71.62471986, 69.06286328, 67.14738432, 60.89947629, 54.80925777, 43.98128998, 29.77145816, 20.08901656, 10.96184517],
  rSTAT5A_rel = [14.72316822, 33.76234229, 36.79985129, 49.71760229, 46.9281201, 47.83657456, 46.92872725, 40.59775294, 43.78366389, 44.45738765, 41.32715926, 41.06273321, 39.23583003, 36.61946054, 34.8937144, 32.21107716]
)

saveat = data.time
tspan = (0.,saveat[end])

u0 = zeros(8)
u0[1] = 207.6*ratio         # STAT5A
u0[3] = 207.6 - 207.6*ratio # STAT5B

ode_jakstat_prob(p) = ODEProblem(stat5_ode, eltype(p).(u0), tspan, p)

# solver algorithm, tolerances
solver_opts = Dict(
    :alg => AutoTsit5(Rodas5P()),
    :reltol => 1e-6,
    :abstol => 1e-8
)

function solve_prob(p)
    _prob = ode_jakstat_prob(p)

    # solution
    sol = solve(_prob, solver_opts[:alg], saveat=saveat, reltol=solver_opts[:reltol],abstol=solver_opts[:abstol]) #save_idxs=[1,2,3,4,5] 
    STAT5A = sol[1,:]
    pApA = sol[2,:]
    STAT5B = sol[3,:]
    pApB = sol[4,:]
    pBpB = sol[5,:]

    # observables
    pSTAT5A_rel = (100 * pApB + 200 * pApA * specC17) ./ (pApB + STAT5A * specC17 + 2 * pApA * specC17)
    pSTAT5B_rel = -(100 * pApB - 200 * pBpB * (specC17 - 1)) ./ ((STAT5B * (specC17 - 1) - pApB) + 2 * pBpB * (specC17 - 1))
    rSTAT5A_rel = (100 * pApB + 100 * STAT5A * specC17 + 200 * pApA * specC17) ./ (2 * pApB + STAT5A * specC17 + 2 * pApA * specC17 - STAT5B * (specC17 - 1) - 2 * pBpB * (specC17 - 1))

    return [pSTAT5A_rel, pSTAT5B_rel, rSTAT5A_rel]
end

p_best = [
    0.026982514033029,      # Epo_degradation_BaF3
    0.0000100067973851508,  # k_exp_hetero
    0.006170228086381,      # k_exp_homo
    0.0163679184468,        # k_imp_hetero
    97749.3794024716,       # k_imp_homo
    15766.5070195731,       # k_phos
    3.85261197844677,       # sd_pSTAT5A_rel
    6.59147818673419,       # sd_pSTAT5B_rel
    3.15271275648527        # sd_rSTAT5A_rel
]

function jakstat_obj(p_init, __p)
  p = exp10.(p_init)

  sim = solve_prob(p)
  σ = p[7:9]
  # loss
  return obj2(sim,data,σ)
end

function obj2(sim,data,σ)
  loss = 0.0
  obs = names(data)[2:end]  

  for i in 1:length(obs)
      loss_i = jakstat_obj_component(sim[i],data[!,i+1],σ[i])
      loss += loss_i
  end
  return loss
end

function jakstat_obj_component(sim,data,σ)
  loss_i = 0.0
  
  for i in eachindex(sim)
          loss_i += ((sim[i]-data[i])/σ)^2 + 2*log(sqrt(2π)*σ)
  end
  return loss_i
end
