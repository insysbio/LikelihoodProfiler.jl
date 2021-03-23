# Overview

**LikelihoodProfiler** is a [Julia](https://julialang.org/downloads/) package for **identifiability analysis** and **confidence intervals** estimation.

## Use cases
Notebooks with use cases can be found in a separate repository: <https://github.com/insysbio/likelihoodprofiler-cases>

 Case | Ref
 ----|----
 PK model with saturation in elimination | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/likelihoodprofiler-cases/master?filepath=notebook/pk_saturation.ipynb)
 Local optim methods comparison | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/likelihoodprofiler-cases/master?filepath=notebook/Derivative-free%20algs%20comparison.ipynb)
 TGF-β Signaling Pathway Model | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/likelihoodprofiler-cases/master?filepath=notebook/TGFb_pathway.ipynb)
 SIR Model. A simple model used as an exercise in identifiability analysis. | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/likelihoodprofiler-cases/master?filepath=notebook/SIR%20Model.ipynb)
 Cancer Taxol Treatment Model  | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/likelihoodprofiler-cases/master?filepath=notebook/taxol_treatment.ipynb)
 
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
  1,             # index of parameter
  f,             # profile function
  :LIN_EXTRAPOL; # method
  loss_crit = 9. # likelihood confidence level
  )
#

# Plot parameter profile x[1]
using Plots
plotly()
plot(res_1)
```

![plot_lin](https://github.com/insysbio/LikelihoodProfiler.jl/blob/master/img/plot_lin.png?raw=true)

## Intro

The reliability and predictability of a **kinetic systems biology (SB) and systems pharmacology (SP) model** depends on the calibration of model parameters. Experimental data can be insufficient to determine all the parameters unambiguously. This results in "non-identifiable" parameters and parameters identifiable within confidence intervals (CIs). The algorithms included in **LikelihoodProfiler** implement Profile Likelihood (PL) [2] method for parameters identification and can be applied to complex SB models. The results of the algorithms can be used to qualify and calibrate parameters or to reduce the model.

## Objective

The package introduces several original algorithms. The default algorithm `:CICO_ONE_PASS` has been developed in accordance with the following principles:

- The algorithms don't require the likelihood function to be differentiable. Hence, derivative-free or global optimization methods can be used to estimate CI endpoints.
- The algorithms are designed to obtain CI endpoints and avoid the calculation of profiles as the most computationally expensive part of the analysis. 
- CI endpoints are estimates with some preset tolerance. Reasonable tolerance setup can also reduce the number of likelihood function calls and speed up the computations. 
- The algorithm are applicable for both finite and infinite CI.

## Methods overview

This algorithms can be applied to complex kinetic models where function differentiability is not guaranteed and each likelihood estimation is computationally expensive.  

The package introduces original "one-pass" algorithm: **Confidence Intervals evaluation by Constrained Optimization** [6]  `:CICO_ONE_PASS` developed by the authors of this package. `:CICO_ONE_PASS` utilizes the **Inequality-based Constrained Optimization** [3-4] for efficient determination of confidence intervals and detection of “non-identifiable” parameters.  

The "multi-pass" methods use extrapolation/interpolation of profile likelihood points: linear (`:LIN_EXTRAPOL`) and quadratic (`:QUADR_EXTRAPOL`) approaches. They are also efficient for both identifiable and non-identifiable parameters.

## References

1. Wikipedia [Identifiability_analysis](https://en.wikipedia.org/wiki/Identifiability_analysis)
2. Kreutz C., et al. Profile Likelihood in Systems Biology. FEBS Journal 280(11), 2564-2571, 2013
3. Steven G. Johnson, The NLopt nonlinear-optimization package, [link](http://ab-initio.mit.edu/nlopt)
4. Andrew R. Conn, Nicholas I. M. Gould, and Philippe L. Toint, "A globally convergent augmented Lagrangian algorithm for optimization with general constraints and simple bounds," SIAM J. Numer. Anal. vol. 28, no. 2, p. 545-572 (1991)
5. Julia: A Fresh Approach to Numerical Computing. Jeff Bezanson, Alan Edelman, Stefan Karpinski and Viral B. Shah (2017) SIAM Review, 59: 65–98
6. Borisov I., Metelkin E. An Algorithm for Practical Identifiability Analysis and Confidence Intervals Evaluation Based on Constrained Optimization. 2018. October. ICSB2018. https://doi.org/10.13140/RG.2.2.18935.06563
