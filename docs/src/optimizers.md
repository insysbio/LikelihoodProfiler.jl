# Note on optimizers

Currently **LikelihoodProfiler** relies on [`NLopt`](https://nlopt.readthedocs.io/en/latest/). 

Internally **LikelihoodProfiler** utilizes `:LN_AUGLAG` algorithm from **NLopt** to construct an augmented objective function with a set of bound constrains. Then the augmented objective function is passed to an optimization algorithm, which is defined by the keyword argument `local_alg` (see [`LikelihoodProfiler.get_interval`](@ref)).

Default `local_alg` = `:LN_NELDERMEAD`, which is the most reliable for the current problem among the derivative-free algorithms. 
Any [`NLopt`](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/) algorithm could be potentially used as `local_alg`.
Nevertheless we performed a series of simulation experiments displaying that some of algorithms may fail at some specific problems.

The following algorithms have shown satisfactory results on test problems (see also `/test/test_grad_algs.jl`):

- `:LN_NELDERMEAD`
- `:LD_MMA`
- `:LD_SLSQP`
- `:LD_CCSAQ`
