
#=
    res1 = get_endpoint(
       [3.0],
       (x) -> (x[1]^2, f_1p(x)),
       :CICO_ONE_PASS,
       :left;
       loss_crit = 8.,
       local_alg=:LD_MMA
    )
=#

algorithms = [
    #derivative-free algorithms
    :LN_NELDERMEAD,
    :LN_SBPLX,
    :LN_COBYLA,
    :LN_BOBYQA,
    :LN_PRAXIS,
    # gradient-based algorithms
    :LD_MMA, # Method of Moving Asymptotes
    :LD_SLSQP, # Sequential Least-Squares Quadratic Programming
    :LD_LBFGS, # Low-storage BFGS
    :LD_TNEWTON_PRECOND_RESTART, # Preconditioned truncated Newton
    :LD_TNEWTON_PRECOND, # Same without restarting
    :LD_TNEWTON_RESTART, # Same without preconditioning
    :LD_TNEWTON, # Same without restarting or preconditioning
    :LD_VAR2, # Shifted limited-memory variable-metric (rank 2)
    :LD_VAR1  # Shifted limited-memory variable-metric (rank 1)
]

res = Dict{Symbol, Union{Nothing,EndPoint}}()

for alg in algorithms
   res[alg] = try get_endpoint(
       [3.0],
       (x) -> (log10(x[1]^2), f_1p(x)),
       :CICO_ONE_PASS,
       :left;
       loss_crit = 8.,
       local_alg=alg,
       max_iter=5000
   )
   catch e
       @warn "$alg error: $e"
   end
end
