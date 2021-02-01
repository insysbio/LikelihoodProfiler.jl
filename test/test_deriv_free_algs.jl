
### Derivative free optimization tests
# The following NLopt derivative-free algorithms are compared
dfo_algorithms = [
    :LN_NELDERMEAD, # Nelder Mead
    :LN_SBPLX, # Subplex (a variant of Nelder-Mead that uses Nelder-Mead on a sequence of subspaces)
    :LN_COBYLA, # Constrained Optimization BY Linear Approximations
    :LN_BOBYQA, # BOBYQA algorithm for bound constrained optimization without derivatives
    :LN_PRAXIS, # "PRAXIS" gradient-free local optimization via the "principal-axis method"
]

@testset "testing dfo algorithms" begin [test_alg(alg; bounds=(-Inf,Inf)) for alg in dfo_algorithms] end
