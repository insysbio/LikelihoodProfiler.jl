# Overview

**LikelihoodProfiler** is a [Julia](https://julialang.org/downloads/) package for **identifiability analysis** and **confidence intervals** evaluation.

## Cases of usage LikelihoodProfiler.jl
 Case | Ref
 ----|----
 PK model with saturation in elimination | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/LikelihoodProfiler.jl/master?filepath=%2Fnotebook%2Fpk_saturation.ipynb)
 Local optim methods comparison | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/LikelihoodProfiler.jl/master?filepath=notebook%2FDerivative-free%20algs%20comparison.ipynb)
 
## Installation


```julia
julia> ]

(v1.2) pkg> add LikelihoodProfiler
```

## Quick start

```julia
using LikelihoodProfiler

# testing profile function
f(x) = 5.0 + (x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2

# Calculate parameters intervals for first parameter component, x[1]
res_1 = get_interval(
  [3., 2., 2.1], # starting point
  1,             # parameter component
  f,             # profile function
  :LIN_EXTRAPOL; # method
  loss_crit = 9. # critical level
  )
#

# Plot parameter profile x[1]
using Plots
plotly()
plot(res_1)
```

![plot_lin](https://github.com/insysbio/LikelihoodProfiler.jl/blob/master/img/plot_lin.png?raw=true)

## Intro

The reliability and predictability of a **kinetic systems biology (SB) and systems pharmacology (SP) model** depends on the calibration of model parameters. Taking into account the lacking of data and the experimental variability the value of any parameter determined unambiguously. This results in characterization of parameter by "confidence intervals" or even "non-identifiable" parameters when the confidence interval is open. The package includes algorithms to perform practical identifiability analysis and evaluation confidence intervals using Profile Likelihood [2] which can be applied to complex SB/SP models. Results of the identifiability analysis can be used to qualify and calibrate parameters or to reduce the model.

## Objective

The package introduces several original algorithms taking into account the following points:

- This algorithm does not assume that the likelihood function is differentiable at any point. This allows using derivation free and global methods of optimization which do not require the calculation of gradients.
- The calculation of likelihood function is the most computationally expensive operation within the others. It becomes critical for large dynamic model used nowadays in systems biology.
- The algorithm should calculate the confidence endpoint with the selected tolerance and must be optimal regarding likelihood function calls. The intermediate (non-endpoint) profile points is not important.
- The algorithm should be stable for calculation both finite and infinite intervals. They should stop immediately (with the corresponding status) if parameter is not identifiable.

## Methods overview

This algorithms can be applied to complex kinetic models where function differentiability is not guaranteed and each likelihood estimation is computationally expensive.  

The package introduces original "one-pass" algorithm: **Confidence Intervals evaluation by Constrained Optimization** [6]  `:CICO_ONE_PASS` developed by the authors of this package. `:CICO_ONE_PASS` utilizes the **Inequality-based Constrained Optimization** [3-4] for efficient determination of confidence intervals and detection of “non-identifiable” parameters.  

The "multi-pass" methods use extrapolation/interpolation of likelihood points to the critical level: linear (`:LIN_EXTRAPOL`) and quadratic (`:QUADR_EXTRAPOL`) approaches. They are also effective for both identifiable and non-identifiable parameters.

## References

1. Wikipedia [Identifiability_analysis](https://en.wikipedia.org/wiki/Identifiability_analysis)
2. Kreutz C., et al. Profile Likelihood in Systems Biology. FEBS Journal 280(11), 2564-2571, 2013
3. Steven G. Johnson, The NLopt nonlinear-optimization package, [link](http://ab-initio.mit.edu/nlopt)
4. Andrew R. Conn, Nicholas I. M. Gould, and Philippe L. Toint, "A globally convergent augmented Lagrangian algorithm for optimization with general constraints and simple bounds," SIAM J. Numer. Anal. vol. 28, no. 2, p. 545-572 (1991)
5. Julia: A Fresh Approach to Numerical Computing. Jeff Bezanson, Alan Edelman, Stefan Karpinski and Viral B. Shah (2017) SIAM Review, 59: 65–98
6. Borisov I., Metelkin E. An Algorithm for Practical Identifiability Analysis and Confidence Intervals Evaluation Based on Constrained Optimization. 2018. October. ICSB2018. https://doi.org/10.13140/RG.2.2.18935.06563
