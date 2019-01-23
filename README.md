# LikelihoodProfiler

[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://insysbio.github.io/LikelihoodProfiler.jl/latest)
[![Travis](https://travis-ci.org/insysbio/LikelihoodProfiler.jl.svg?branch=master)](https://travis-ci.org/insysbio/LikelihoodProfiler.jl)
[![AppVeyor](https://ci.appveyor.com/api/projects/status/8qdhx23slm9qc0m2?svg=true)](https://ci.appveyor.com/project/ivborissov/likelihoodprofiler-jl)
[![GitHub release](https://img.shields.io/github/release/insysbio/LikelihoodProfiler.jl.svg)](https://github.com/insysbio/LikelihoodProfiler.jl/releases/)
[![GitHub license](https://img.shields.io/github/license/insysbio/LikelihoodProfiler.jl.svg)](https://github.com/insysbio/LikelihoodProfiler.jl/blob/master/LICENSE)
[![DOI:10.13140/RG.2.2.18935.06563](https://zenodo.org/badge/DOI/10.13140/RG.2.2.18935.06563.svg)](https://doi.org/10.13140/RG.2.2.18935.06563)

**LikelihoodProfiler** is a Julia package for **identifiability analysis** and confidence intervals evaluation.

# Installation

```
julia> import Pkg   # if you are on Julia 0.7

julia> Pkg.add(PackageSpec(url="https://github.com/insysbio/LikelihoodProfiler.jl.git"))

julia> using LikelihoodProfiler
```

# Quick start

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
        :CICO_ONE_PASS;
        loss_crit = 9.
    ) for i in 1:3]

# Plot parameter profile x[1]
using Plots
plotly()
plot(res[1])
```
