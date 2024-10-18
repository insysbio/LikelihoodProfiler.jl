all_algorithms_loss = [
    # good
    (algorithm = :LN_NELDERMEAD, skip = [:f_1p, :f_2p_1im]), # Nelder Mead
    # unstable
    (algorithm = :LN_PRAXIS, skip = [:f_1p]), # "PRAXIS" gradient-free local optimization via the "principal-axis method"
    (algorithm = :LN_SBPLX, skip = []), # Subplex (a variant of Nelder-Mead that uses Nelder-Mead on a sequence of subspaces)
    # errors
    (algorithm = :LN_NEWUOA, skip = [:f_1p,:f_4p_2im,:f_3p_im,:f_3p_1im,:f_3p_1im_dep,:f_1p_ex,:f_4p_3im,:f_2p_1im,:f_2p]),
    (algorithm = :LN_BOBYQA, skip = [:f_1p,:f_3p_1im,:f_4p_2im,:f_1p_ex,:f_4p_3im,:f_5p_3im,:f_2p_1im,:f_2p]), # BOBYQA algorithm for bound constrained optimization without derivatives
    (algorithm = :LN_COBYLA, skip = []), # Constrained Optimization BY Linear Approximations
]
all_algorithms_loss_log = [
    # good
    (algorithm = :LN_NELDERMEAD, skip = [:f_1p,:f_2p_1im]), # Nelder Mead
    # unstable
    (algorithm = :LN_PRAXIS, skip = [:f_1p]), # "PRAXIS" gradient-free local optimization via the "principal-axis method"
    (algorithm = :LN_SBPLX, skip = []), # Subplex (a variant of Nelder-Mead that uses Nelder-Mead on a sequence of subspaces)
    # errors
    (algorithm = :LN_NEWUOA, skip = [:f_1p,:f_4p_2im,:f_3p_im]), # Unconstrained Optimization BY Quadratic Approximations
    (algorithm = :LN_BOBYQA, skip = [:f_1p_ex,:f_3p_im,:f_5p_3im]), # BOBYQA algorithm for bound constrained optimization without derivatives
    (algorithm = :LN_COBYLA, skip = [:f_1p_ex,:f_3p_im]), # Constrained Optimization BY Linear Approximations
]
all_algorithms_scan = [
    # good
    (algorithm = :LN_NELDERMEAD, skip = []), # Nelder Mead
    # unstable
    (algorithm = :LN_PRAXIS, skip = [:f_1p]), # "PRAXIS" gradient-free local optimization via the "principal-axis method"
    (algorithm = :LN_SBPLX, skip = []), # Subplex (a variant of Nelder-Mead that uses Nelder-Mead on a sequence of subspaces)
    # errors
    (algorithm = :LN_NEWUOA, skip = [:f_1p,:f_4p_2im,:f_3p_im,:f_3p_1im_dep,:f_1p_ex,:f_2p_1im,:f_2p,:f_5p_3im]),
    (algorithm = :LN_BOBYQA, skip = [:f_3p_im]), # BOBYQA algorithm for bound constrained optimization without derivatives
    (algorithm = :LN_COBYLA, skip = []), # Constrained Optimization BY Linear Approximations
]
all_algorithms_scan_log = [
    # good
    (algorithm = :LN_NELDERMEAD, skip = []), # Nelder Mead
    # unstable
    (algorithm = :LN_PRAXIS, skip = [:f_1p]), # "PRAXIS" gradient-free local optimization via the "principal-axis method"
    (algorithm = :LN_SBPLX, skip = []), # Subplex (a variant of Nelder-Mead that uses Nelder-Mead on a sequence of subspaces)
    # errors
    (algorithm = :LN_NEWUOA, skip = [:f_1p,:f_3p_1im,:f_3p_1im_dep,:f_4p_2im,:f_1p_ex,:f_2p_1im,:f_2p,:f_3p_im]),
    (algorithm = :LN_BOBYQA, skip = [:f_3p_1im,:f_3p_1im_dep,:f_3p_im,:f_4p_3im,:f_5p_3im]), # BOBYQA algorithm for bound constrained optimization without derivatives
    (algorithm = :LN_COBYLA, skip = []), # Constrained Optimization BY Linear Approximations
]

@testset "loss" begin
    [test_alg_optimal(alg; loss_tol = 1e-3) for alg in all_algorithms_loss] 
end

@testset ":log" begin
    [test_alg_optimal(alg; loss_tol = 1e-3, scale = :log, bounds = (0.,Inf)) for alg in all_algorithms_loss_log] 
end

@testset "scan" begin
    [test_alg_optimal(alg; scan_tol = 1e-4) for alg in all_algorithms_scan] 
end

@testset "scan :log" begin
    [test_alg_optimal(alg; scan_tol = 1e-4, scale = :log, bounds = (0.,Inf)) for alg in all_algorithms_scan_log] 
end
