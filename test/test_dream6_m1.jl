include("./dream6_m1/dream6_m1_@macro.jl")

# using JSON

using NLopt, LikelihoodProfiler

# Scale types
mixed_scale = fill(true, length(p0))
mixed_scale[[15,17,19,21,23,25,27,29]] = false

# Bounds
bounds_params = fill([1e-12, 1e12], length(p0))
bounds_params[[15,17,19,21,23,25,27,29]] = fill([0.5, 9.], 8)

# fit optimum
opt_0 = Opt(
    :LN_NELDERMEAD, # :LN_SBPLX
    length(p0)
)

min_objective!(opt_0, obj)
ftol_abs!(opt_0, 1e-3)
# bound constraints
lb = minimum.(bounds_params)
ub = maximum.(bounds_params)
lower_bounds!(opt_0, lb)
upper_bounds!(opt_0, ub)
(obj_0, p_0, ret_0) = optimize(opt_0, p0)

# calculate critical level
obj_crit = obj_0 + chisqinvcdf(1, 0.95)

# scan bounds
scan_bounds = fill([1e-6, 1e6], length(p0))
scan_bounds[[15,17,19,21,23,25,27,29]] = fill([0.6, 8.], 8)

# calculate intervals for all
max_iter = 100000
result = []
for i = 1:length(p0)
    println()
    println("START # $i")
    out_inq = params_intervals(
        p_0,
        i,
        obj_crit,
        obj,
        logscale = mixed_scale,
        bounds = bounds_params,
        scan_bound = scan_bounds[i],
        max_iter = max_iter
    )
    println(
        "FINISH #$i, iteration count: ",
        round(out_inq.count_evals[1]/max_iter*100, 1), "%   ",
        round(out_inq.count_evals[2]/max_iter*100, 1), "%"
    )
    push!(result, out_inq)
end

result

open("result.json", "w") do f
    JSON.print(f, result, 4)
end

open("result.csv", "w") do f
    CSV.print(f, result)
end

[println(r.loss_final) for r in result]
