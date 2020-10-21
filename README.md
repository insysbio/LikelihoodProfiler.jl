# LikelihoodProfiler

[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://insysbio.github.io/LikelihoodProfiler.jl/latest)
[![Travis](https://travis-ci.org/insysbio/LikelihoodProfiler.jl.svg?branch=master)](https://travis-ci.org/insysbio/LikelihoodProfiler.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/ntk7f1lpjct58n6t/branch/master?svg=true)](https://ci.appveyor.com/project/metelkin/likelihoodprofiler-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/insysbio/LikelihoodProfiler.jl/badge.svg?branch=master)](https://coveralls.io/github/insysbio/LikelihoodProfiler.jl?branch=master)
[![GitHub release](https://img.shields.io/github/release/insysbio/LikelihoodProfiler.jl.svg)](https://github.com/insysbio/LikelihoodProfiler.jl/releases/)
[![GitHub license](https://img.shields.io/github/license/insysbio/LikelihoodProfiler.jl.svg)](https://github.com/insysbio/LikelihoodProfiler.jl/blob/master/LICENSE)
[![DOI:10.13140/RG.2.2.10306.94409](https://zenodo.org/badge/DOI/10.13140/RG.2.2.10306.94409.svg)](https://doi.org/10.13140/RG.2.2.10306.94409)

**LikelihoodProfiler** is a [Julia language](https://julialang.org/downloads/) package for [identifiability analysis](https://en.wikipedia.org/wiki/Identifiability_analysis) and confidence intervals evaluation.

See [documentation](https://insysbio.github.io/LikelihoodProfiler.jl/latest/).

## Cases
Cases notebooks have been removed to separate repository: <https://github.com/insysbio/likelihoodprofiler-cases>

 Case | Ref
 ----|----
 PK model with saturation in elimination | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/likelihoodprofiler-cases/master?filepath=notebook/pk_saturation.ipynb)
 Local optim methods comparison | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/likelihoodprofiler-cases/master?filepath=notebook/Derivative-free%20algs%20comparison.ipynb)
 TGF-β Signaling Pathway Model | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/likelihoodprofiler-cases/master?filepath=notebook/TGFb_pathway.ipynb)
 SIR Model. A simple model used as an exercise in identifiability analysis. | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/likelihoodprofiler-cases/master?filepath=notebook/SIR%20Model.ipynb)
 Cancer Taxol Treatment Model  | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/likelihoodprofiler-cases/master?filepath=notebook/taxol_treatment.ipynb)

# Installation

```julia
julia> ]

(v1.2) pkg> add LikelihoodProfiler
```

# Quick start

```julia
using LikelihoodProfiler

# testing profile function
f(x) = 5.0 + (x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2

# Calculate parameters intervals for first parameter component, x[1]
res_1 = get_interval(
  [3., 2., 2.1], # starting point
  1,             # parameter component to analyze
  f,             # profile function
  :LIN_EXTRAPOL; # method
  loss_crit = 9. # critical level of loss function
  )
#

# Plot parameter profile x[1]
using Plots
plotly()
plot(res_1)
```

![Plot Linear](img/plot_lin.png?raw=true)
