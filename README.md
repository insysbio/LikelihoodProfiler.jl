# LikelihoodProfiler

[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://insysbio.github.io/LikelihoodProfiler.jl/latest)
[![Build status](https://github.com/insysbio/LikelihoodProfiler.jl/workflows/CI/badge.svg)](https://github.com/insysbio/LikelihoodProfiler.jl/actions?query=workflow%3ACI)
[![Coverage Status](https://coveralls.io/repos/github/insysbio/LikelihoodProfiler.jl/badge.svg?branch=master)](https://coveralls.io/github/insysbio/LikelihoodProfiler.jl?branch=master)
[![GitHub release](https://img.shields.io/github/release/insysbio/LikelihoodProfiler.jl.svg)](https://github.com/insysbio/LikelihoodProfiler.jl/releases/)
[![GitHub license](https://img.shields.io/github/license/insysbio/LikelihoodProfiler.jl.svg)](https://github.com/insysbio/LikelihoodProfiler.jl/blob/master/LICENSE)
[![DOI:10.1371/journal.pcbi.1008495](https://zenodo.org/badge/DOI/10.1371/journal.pcbi.1008495.svg)](https://doi.org/10.1371/journal.pcbi.1008495)

**LikelihoodProfiler** is a [Julia language](https://julialang.org/downloads/) package for [identifiability analysis](https://en.wikipedia.org/wiki/Identifiability_analysis) and confidence intervals estimation.

See [documentation](https://insysbio.github.io/LikelihoodProfiler.jl/latest/).

## Use cases
Notebooks with use cases can be found in a separate repository: <https://github.com/insysbio/likelihoodprofiler-cases>

 Case | Ref
 ----|----
 PK model with saturation in elimination | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/likelihoodprofiler-cases/master?filepath=notebook/pk_saturation.ipynb)
 Local optim methods comparison | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/likelihoodprofiler-cases/master?filepath=notebook/Derivative-free%20algs%20comparison.ipynb)
 TGF-β Signaling Pathway Model | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/likelihoodprofiler-cases/master?filepath=notebook/TGFb_pathway.ipynb)
 SIR Model. A simple model used as an exercise in identifiability analysis. | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/likelihoodprofiler-cases/master?filepath=notebook/SIR%20Model.ipynb)
 Cancer Taxol Treatment Model  | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/likelihoodprofiler-cases/master?filepath=notebook/taxol_treatment.ipynb)
 STAT5 Dimerization Model  | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/insysbio/likelihoodprofiler-cases/master?filepath=notebook/STAT5%20Dimerization.ipynb)
 

## Installation

```julia
julia> ]

(v1.7) pkg> add LikelihoodProfiler
```

## Quick start

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

## License

[MIT](./LICENSE) Public license

## How to cite

**Borisov I, Metelkin E** (2020) *Confidence intervals by constrained optimization—An algorithm and software package for practical identifiability analysis in systems biology.* PLoS Comput Biol 16(12): e1008495.

Ref: <https://doi.org/10.1371/journal.pcbi.1008495>
