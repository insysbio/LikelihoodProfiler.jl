# Overview

**LikelihoodProfiler** is a Julia package for **identifiability analysis** and confidence intervals evaluation.

## Installation

**Julia** [download page](https://julialang.org/downloads/).

Currently supported **Julia** versions are 0.7, 1.0

To install the package from `REPL`
```
julia> import Pkg   # if you are on Julia 0.7, 1.0

julia> Pkg.add(PackageSpec(url="https://github.com/insysbio/LikelihoodProfiler.jl.git"))

julia> using LikelihoodProfiler
```

## Quick start

```
using LikelihoodProfiler

# Likelihood function
f(x) = 5.0 + (x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2

# Calculate parameters intervals for x[1], x[2], x[3]
res = [
    get_interval(
        [3., 2., 2.1],
        i,
        f,
        :LIN_EXTRAPOL;
        loss_crit = 9.
    ) for i in 1:3]

# Plot parameter profile x[1]
using Plots
plotly()
plot(res[1])
```

![](https://github.com/insysbio/LikelihoodProfiler.jl/blob/master/img/plot_lin.png?raw=true)

## Objective

The reliability and predictability of a **kinetic systems biology (SB) model** depends on the calibration of model parameters. Experimental data can be insufficient to determine all the parameters unambiguously. This results in “non-identifiable” parameters and parameters identifiable within confidence intervals. The package includes algorithms for parameters identification using Profile Likelihood [1] method which can be applied to complex SB models. Results of the identifiability analysis can be used to qualify and calibrate parameters or to reduce the model.


## Package Overview

This packages provides a number of algorithms for [identifiability analysis](https://en.wikipedia.org/wiki/Identifiability_analysis) and confidence intervals evaluation based on **Profile Likelihood** approach. Along with linear extrapolation (`:LIN_EXTRAPOL`) and quadratic extrapolation (`:QUADR_EXTRAPOL`) the package introduces **Confidence Intervals evaluation by Constrained Optimization** algorithm (`:CICO_ONE_PASS`) developed by the authors of this package.
`:CICO_ONE_PASS` utilizes the **Inequality-based Constrained Optimization** [2, 3] for efficient determination of confidence intervals and detection of “non-identifiable” parameters. This algorithm does not assume that the likelihood function is differentiable or can be calculated for any given parameters set. This algorithm can be applied to complex kinetic models where function differentiability is not guaranteed and each likelihood estimation is computationally expensive.  


## References

1. Kreutz C., et al. Profile Likelihood in Systems Biology. FEBS Journal 280(11), 2564-2571, 2013
2. Steven G. Johnson, The NLopt nonlinear-optimization package, [link](http://ab-initio.mit.edu/nlopt)
3. Andrew R. Conn, Nicholas I. M. Gould, and Philippe L. Toint, "A globally convergent augmented Lagrangian algorithm for optimization with general constraints and simple bounds," SIAM J. Numer. Anal. vol. 28, no. 2, p. 545-572 (1991)
4. Julia: A Fresh Approach to Numerical Computing. Jeff Bezanson, Alan Edelman, Stefan Karpinski and Viral B. Shah (2017) SIAM Review, 59: 65–98
