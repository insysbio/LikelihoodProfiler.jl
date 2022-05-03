
### Gradient-based optimization tests
# The following NLopt gradient -based algorithms are compared

grad_algorithms = [
    # good
    #:LD_MMA, # Method of Moving Asymptotes #tmp removed from tests
    (algorithm = :LD_SLSQP, skip = []), # Sequential Least-Squares Quadratic Programming
    (algorithm = :LD_CCSAQ, skip = []), # Conservative convex separable approximation
    # errors
    #(algorithm = :LD_LBFGS, skip = []), # Low-storage BFGS
    #(algorithm = :LD_TNEWTON_PRECOND_RESTART, skip = []), # Preconditioned truncated Newton
    #(algorithm = :LD_TNEWTON_PRECOND, skip = []), # Same without restarting
    #(algorithm = :LD_TNEWTON_RESTART, skip = []), # Same without preconditioning
    #(algorithm = :LD_TNEWTON, skip = []), # Same without restarting or preconditioning
    #(algorithm = :LD_VAR2, skip = []), # Shifted limited-memory variable-metric (rank 2)
    #(algorithm = :LD_VAR1, skip = [])  # Shifted limited-memory variable-metric (rank 1)
]

[test_alg(alg; bounds=(-1e10,1e10)) for alg in grad_algorithms]
