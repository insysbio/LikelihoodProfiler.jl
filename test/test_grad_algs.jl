
### Gradient-based optimization tests
# The following NLopt gradient -based algorithms are compared
grad_algorithms = [
    #:LD_MMA, # Method of Moving Asymptotes #tmp removed from tests
    :LD_SLSQP, # Sequential Least-Squares Quadratic Programming
    :LD_CCSAQ, # Conservative convex separable approximation
    # errors
    #:LD_LBFGS, # Low-storage BFGS
    #:LD_TNEWTON_PRECOND_RESTART, # Preconditioned truncated Newton
    #:LD_TNEWTON_PRECOND, # Same without restarting
    #:LD_TNEWTON_RESTART, # Same without preconditioning
    #:LD_TNEWTON, # Same without restarting or preconditioning
    #:LD_VAR2, # Shifted limited-memory variable-metric (rank 2)
    #:LD_VAR1  # Shifted limited-memory variable-metric (rank 1)
]

[test_alg(alg; bounds=(-1e10,1e10)) for alg in grad_algorithms]
