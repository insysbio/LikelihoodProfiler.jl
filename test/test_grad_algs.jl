
### Gradient-based optimization tests
# The following NLopt gradient -based algorithms are compared

grad_algorithms_autodiff = [
    # good
    #:LD_MMA, # Method of Moving Asymptotes #tmp removed from tests
    (algorithm = :LD_SLSQP, skip = []), # Sequential Least-Squares Quadratic Programming
    (algorithm = :LD_CCSAQ, skip = []), # Conservative convex separable approximation
    # errors
    (algorithm = :LD_LBFGS, skip = [:f_3p_1im, :f_4p_2im, :f_4p_3im, :f_5p_3im, :f_2p]), # Low-storage BFGS
    (algorithm = :LD_TNEWTON_PRECOND_RESTART, skip = [:f_3p_1im, :f_4p_2im, :f_4p_3im, :f_5p_3im, :f_2p]), # Preconditioned truncated Newton
    (algorithm = :LD_TNEWTON_PRECOND, skip = [:f_3p_1im, :f_4p_3im, :f_5p_3im]), # Same without restarting
    (algorithm = :LD_TNEWTON_RESTART, skip = [:f_3p_1im, :f_4p_2im, :f_4p_3im, :f_2p]), # Same without preconditioning
    (algorithm = :LD_TNEWTON, skip = [:f_3p_1im, :f_4p_3im]), # Same without restarting or preconditioning
    (algorithm = :LD_VAR2, skip = [:f_3p_1im, :f_4p_2im, :f_4p_3im, :f_5p_3im, :f_2p]), # Shifted limited-memory variable-metric (rank 2)
    (algorithm = :LD_VAR1, skip = [:f_3p_1im, :f_4p_2im, :f_4p_3im, :f_5p_3im, :f_2p])  # Shifted limited-memory variable-metric (rank 1)
]

[test_alg(alg; bounds=(-1e10,1e10), loss_grad=:AUTODIFF) for alg in grad_algorithms_autodiff]

grad_algorithms_finite = [
    # good
    #:LD_MMA, # Method of Moving Asymptotes #tmp removed from tests
    (algorithm = :LD_SLSQP, skip = []), # Sequential Least-Squares Quadratic Programming
    (algorithm = :LD_CCSAQ, skip = [:f_5p_3im]), # Conservative convex separable approximation
    # errors
    (algorithm = :LD_LBFGS, skip = [:f_3p_1im, :f_4p_3im, :f_5p_3im]), # Low-storage BFGS
    (algorithm = :LD_TNEWTON_PRECOND_RESTART, skip = [:f_3p_1im, :f_3p_1im_dep, :f_4p_3im, :f_5p_3im]), # Preconditioned truncated Newton
    (algorithm = :LD_TNEWTON_PRECOND, skip = [:f_3p_1im, :f_4p_3im, :f_5p_3im]), # Same without restarting
    (algorithm = :LD_TNEWTON_RESTART, skip = [:f_3p_1im, :f_4p_3im, :f_5p_3im]), # Same without preconditioning
    (algorithm = :LD_TNEWTON, skip = [:f_3p_1im, :f_4p_3im, :f_5p_3im]), # Same without restarting or preconditioning
    (algorithm = :LD_VAR2, skip = [:f_3p_1im, :f_4p_2im, :f_4p_3im, :f_5p_3im, :f_2p]), # Shifted limited-memory variable-metric (rank 2)
    (algorithm = :LD_VAR1, skip = [:f_3p_1im, :f_4p_2im, :f_4p_3im, :f_5p_3im, :f_2p])  # Shifted limited-memory variable-metric (rank 1)
]

[test_alg(alg; bounds=(-1e10,1e10), loss_grad=:FINITE) for alg in grad_algorithms_finite]