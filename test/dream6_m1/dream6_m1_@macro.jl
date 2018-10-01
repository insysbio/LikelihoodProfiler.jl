# Pkg.add("DiffEqMonteCarlo")
# Pkg.add("DiffEqParamEstim")
# Pkg.add("CSV")
# Pkg.add("Sundials")
# Pkg.add("IterableTables")
# Pkg.add("DifferentialEquations")

using DiffEqBase, OrdinaryDiffEq, ParameterizedFunctions,
    DiffEqParamEstim, CSV, StatsFuns, DataFrames, LSODA,
    IterableTables, Distributions, DataStructures

# not used
# using  DiffEqMonteCarlo, Sundials, DifferentialEquations

#using PrintLog
#@printlog "DREAM6_identification.log"

# Input data
exp_data = CSV.read("./test/dream6_m1/exp_data.csv", types=Dict("t"=>Float64), nullable=false)
u_data = CSV.read("./test/dream6_m1/u_data.csv", types=Dict("values"=>Float64), nullable=false)
p_data = CSV.read("./test/dream6_m1/p_data.csv", types=Dict("values"=>Float64), nullable=false)

# Variables and parameters names and data
u_names = OrderedDict(Symbol(u_data[:names][i]) => i for i in eachindex(u_data[:names]))
p_names = OrderedDict(Symbol(p_data[:names][i]) => i for i in eachindex(p_data[:names]))
u0 = u_data[:values]
p0 = p_data[:values]

# Constants
const pp1_mrna_degradation_rate = 1.0
const pp2_mrna_degradation_rate = 1.0
const pp3_mrna_degradation_rate = 1.0
const pp4_mrna_degradation_rate = 1.0
const pp5_mrna_degradation_rate = 1.0
const pp6_mrna_degradation_rate = 1.0

# Power function
pow(x, y) =  x >= 0. ? x^y : 0.0

# DREAM6 Model1
Model1 = @ode_def_bare DREAM6 begin
    rs1 = 1. / (1. + pow(p6 / v5_Kd, v5_h))
    rs2 = 1. / (1. + pow(p5 / v8_Kd, v8_h))
    rs3 = 1. / (1. + pow(p4 / v6_Kd, v6_h))
    rs4 = 1. / (1. + pow(p2 / v4_Kd, v4_h))
    rs5 = 1. / (1. + pow(p4 / v7_Kd, v7_h))

    as1 = pow(p1 / v2_Kd, v2_h) / (1. + pow(p1 / v2_Kd, v2_h))
    as2 = pow(p1 / v1_Kd, v1_h) / (1. + pow(p1 / v1_Kd, v1_h))
    as3 = pow(p1 / v3_Kd, v3_h) / (1. + pow(p1 / v3_Kd, v3_h))

    cod1 = pro1_strength
    cod2 = pro2_strength * as1 * rs1
    cod3 = pro3_strength * as3 * rs4
    cod4 = pro4_strength * as2 * rs2
    cod5 = pro5_strength * rs3
    cod6 = pro6_strength * rs5

    # V = zeros(24)
    V_1 = cod1
    V_2 = pp1_mrna_degradation_rate * pp1_mrna
    V_3 = rbs1_strength * pp1_mrna
    V_4 = p_degradation_rate * p1
    V_5 = cod5
    V_6 = pp5_mrna_degradation_rate * pp5_mrna
    V_7 = rbs5_strength * pp5_mrna
    V_8 = p_degradation_rate * p5
    V_9 = cod2
    V_10 = pp2_mrna_degradation_rate * pp2_mrna
    V_11 = rbs2_strength * pp2_mrna
    V_12 = p_degradation_rate * p2
    V_13 = cod6
    V_14 = pp6_mrna_degradation_rate * pp6_mrna
    V_15 = rbs6_strength * pp6_mrna
    V_16 = p_degradation_rate * p6
    V_17 = cod3
    V_18 = pp3_mrna_degradation_rate * pp3_mrna
    V_19 = rbs3_strength * pp3_mrna
    V_20 = p_degradation_rate * p3
    V_21 = cod4
    V_22 = pp4_mrna_degradation_rate * pp4_mrna
    V_23 = rbs4_strength * pp4_mrna
    V_24 = p_degradation_rate * p4

    # for pp1_mrna
    dpp1_mrna = V_1 - V_2
    # for p1
    dp1 = V_3 - V_4
    # for pp5_mrna
    dpp5_mrna = V_5 - V_6
    # for p5
    dp5 = V_7 - V_8
    # for pp2_mrna
    dpp2_mrna = V_9 - V_10
    # for p2
    dp2 = V_11 - V_12
    # for pp6_mrna
    dpp6_mrna = V_13 - V_14
    # for p6
    dp6 = V_15 - V_16
    # for pp3_mrna
    dpp3_mrna = V_17 - V_18
    # for p3
    dp3 = V_19 - V_20
    # for pp4_mrna
    dpp4_mrna = V_21 - V_22
    # for p4
    dp4 = V_23 - V_24
end p_degradation_rate pro1_strength pro2_strength pro3_strength pro4_strength pro5_strength pro6_strength rbs1_strength rbs2_strength rbs3_strength rbs4_strength rbs5_strength rbs6_strength v1_Kd v1_h v2_Kd v2_h v3_Kd v3_h v4_Kd v4_h v5_Kd v5_h v6_Kd v6_h v7_Kd v7_h v8_Kd v8_h

#prob = ODEProblem(Model1, u0, (0., 30.), params)

#sol = solve(ODEProblem(Model1, u0, (0., 20.), params), Tsit5(), saveat=exp_data[:t])

# Loss function
function ofv_loss(sol)
    loss = 0.0

    if sol.retcode != :Success #any((s.retcode != :Success for s in sol))
        loss = Inf
        throw(ForcedStop())
        println(sol.retcode, sol)
    else
        # println(sol.prob.p)
        for u in names(exp_data)[2:end]
            loss += sum(((sol(exp_data[:t], idxs = u_names[u]).u - exp_data[u]).^2) ./ ((0.1)^2+(0.2*sol(exp_data[:t],idxs = u_names[u]).u).^2))
        end
    end
    return loss
end

# Obj function to optimize
obj = build_loss_objective(
    ODEProblem(Model1, u0, (0., 20.), p0),
    lsoda(), # Rosenbrock23() Tsit5()
    ofv_loss,
    verbose = false,
    maxiters = 1e6,
    saveat = exp_data[:t]
)

#=
using PrintLog
@printlog "DREAM6_param_SBPLX.log"

for scale in (direct_scale, log_scale, mixed_scale), bounds in (bounds_params_open, bounds_params_closed, bounds_params_half_open), alg in (:LN_SBPLX, :LN_COBYLA)

    scale == direct_scale && println("direct_scale")
    scale == log_scale && println("log_scale")
    scale == mixed_scale && println("mixed_scale")

    bounds == bounds_params_open && println("bounds_params_open")
    bounds == bounds_params_closed && println("bounds_params_closed")
    bounds == bounds_params_half_open && println("bounds_params_half_open")
    println(alg)

    params_intervals(p0, 1, obj0, obj; logscale = scale, local_alg = alg, bounds_id = bounds_id[1], bounds_params = bounds)
end

for i in 1:length(p_names)
    params_intervals(p0, i, obj0, obj; logscale = direct_scale, local_alg = :LN_SBPLX, bounds_id = bounds_id[i], bounds_params = bounds_params_half_open)
end
=#
