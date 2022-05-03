### Derivative free optimization tests

using LikelihoodProfiler

# The following NLopt derivative-free algorithms are compared
dfo_algorithms = [
    # good
    (algorithm = :LN_NELDERMEAD, skip = []), # Nelder Mead
    # unstable
    (algorithm = :LN_PRAXIS, skip = [:f_3p_1im, :f_4p_3im]), # "PRAXIS" gradient-free local optimization via the "principal-axis method"
    # errors
    (algorithm = :LN_SBPLX, skip = [:f_3p_1im, :f_3p_1im_dep]), # Subplex (a variant of Nelder-Mead that uses Nelder-Mead on a sequence of subspaces)
    (algorithm = :LN_NEWUOA, skip = [:f_1p, :f_3p_1im, :f_4p_3im, :f_5p_3im, :f_2p_1im]),
    (algorithm = :LN_BOBYQA, skip = [:f_1p, :f_3p_1im, :f_1p_ex, :f_4p_3im, :f_5p_3im, :f_2p_1im]), # BOBYQA algorithm for bound constrained optimization without derivatives
    #(algorithm = :LN_COBYLA, skip = []), # Constrained Optimization BY Linear Approximations
]

[test_alg(alg; bounds=(-Inf,Inf)) for alg in dfo_algorithms]
