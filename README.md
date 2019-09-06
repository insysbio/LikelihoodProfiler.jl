# LikelihoodProfiler

[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://insysbio.github.io/LikelihoodProfiler.jl/latest)
[![Travis](https://travis-ci.org/insysbio/LikelihoodProfiler.jl.svg?branch=master)](https://travis-ci.org/insysbio/LikelihoodProfiler.jl)
[![AppVeyor](https://ci.appveyor.com/api/projects/status/github/insysbio/LikelihoodProfiler.jl?branch=master&svg=true)](https://ci.appveyor.com/project/metelkin/likelihoodprofiler-jl)
[![Coverage Status](https://coveralls.io/repos/github/insysbio/LikelihoodProfiler.jl/badge.svg?branch=master)](https://coveralls.io/github/insysbio/LikelihoodProfiler.jl?branch=master)
[![GitHub release](https://img.shields.io/github/release/insysbio/LikelihoodProfiler.jl.svg)](https://github.com/insysbio/LikelihoodProfiler.jl/releases/)
[![GitHub license](https://img.shields.io/github/license/insysbio/LikelihoodProfiler.jl.svg)](https://github.com/insysbio/LikelihoodProfiler.jl/blob/master/LICENSE)
[![DOI:10.13140/RG.2.2.18935.06563](https://zenodo.org/badge/DOI/10.13140/RG.2.2.18935.06563.svg)](https://doi.org/10.13140/RG.2.2.18935.06563)

**LikelihoodProfiler** is a [Julia language](https://julialang.org/downloads/) package for [identifiability analysis](https://en.wikipedia.org/wiki/Identifiability_analysis) and confidence intervals evaluation.

See [documentation](https://insysbio.github.io/LikelihoodProfiler.jl/latest/).

# Installation

Currently supported **Julia** versions are 0.7, 1.0

```
julia> import Pkg   # if you are on Julia 0.7, 1.0

julia> Pkg.add("LikelihoodProfiler")
```

# Quick start

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

![Plot Linear](img/plot_lin.png?raw=true)
