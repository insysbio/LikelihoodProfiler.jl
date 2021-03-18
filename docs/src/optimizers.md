# Note on optimizers

Currently **LikelihoodProfiler** relies on [`NLopt`](https://nlopt.readthedocs.io/en/latest/). 

Internally **LikelihoodProfiler** utilizes `:LN_AUGLAG` algorithm from **NLopt** to construct an augmented objective function. Then the augmented objective function with no constraints is passed to an optimization algorithm, which is defined by the keyword argument `local_alg` (see [`LikelihoodProfiler.get_interval`](@ref)). Default is `:LN_NELDERMEAD`, which is the most reliable for the current problem among the derivative-free algorithms. The following gradien based algorithms have also shown good relsults on test problems (see /test/test_grad_algs.jl): `:LD_MMA`, `:LD_SLSQP`, `:LD_CCSAQ`.