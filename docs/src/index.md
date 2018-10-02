# Overview

**LikelihoodProfiler** is a Julia package for **identifiability analysis** and confidence intervals evaluation.

## Installation
**Julia** [download page](https://julialang.org/downloads/).

Currently supported **Julia** versions are 0.6, 0.7

To install the package from `REPL`
```
julia> import Pkg   # if you are on Julia 0.7

julia> Pkg.add(PackageSpec(url="https://github.com/insysbio/LikelihoodProfiler.jl.git"))

julia> using LikelihoodProfiler
```

## Objective

The reliability and predictability of a **kinetic systems biology (SB) model** depends on the calibration of model parameters. Experimental data can be insufficient to determine all the parameters unambiguously. This results in “non-identifiable” parameters and parameters identifiable within confidence intervals. The proposed algorithm is a practical implementation of Profile Likelihood [1] method for parameters identification which can be applied to complex SB models. The results of this algorithm can be used to qualify and calibrate parameters or to reduce the model.

## Algorithm

The proposed algorithm for **Profile Likelihood** method addresses the disadvantages and restrictions of the root-finding algorithms with regard to the above problem and utilizes the **Inequality-based Constrained Optimization** [2, 3] for efficient determination of confidence intervals and detection of “non-identifiable” parameters. This algorithm does not assume that the likelihood function is differentiable or can be calculated for any given parameters set. This algorithm can be applied to complex kinetic models where function differentiability is not guaranteed and each likelihood estimation is computationally expensive.  The algorithm was tested for the set of kinetic models and it is distributed as a software package based on Julia Programming Language [4]. The package includes tools for **parameters identifiability analysis**, **confidence intervals evaluation** and **results visualization**.


## References

1. Kreutz C., et al. Profile Likelihood in Systems Biology. FEBS Journal 280(11), 2564-2571, 2013
2. Steven G. Johnson, The NLopt nonlinear-optimization package, [link](http://ab-initio.mit.edu/nlopt)
3. Andrew R. Conn, Nicholas I. M. Gould, and Philippe L. Toint, "A globally convergent augmented Lagrangian algorithm for optimization with general constraints and simple bounds," SIAM J. Numer. Anal. vol. 28, no. 2, p. 545-572 (1991)
4. Julia: A Fresh Approach to Numerical Computing. Jeff Bezanson, Alan Edelman, Stefan Karpinski and Viral B. Shah (2017) SIAM Review, 59: 65–98
