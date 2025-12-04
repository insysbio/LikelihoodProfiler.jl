## Overview

`LikelihoodProfiler.jl` is a [Julia](https://julialang.org/downloads/) package for practical identifiability analysis and confidence intervals estimation using the profile likelihood approach. The package provides a unified interface for various profile likelihood methods, including optimization-based **OptimizationProfiler** and integration-based profiles **IntegrationProfiler**, CI endpoints search **CICOProfiler**, and more.

## Who is this package for?

`LikelihoodProfiler.jl` is intended for researchers and practitioners working with **maximum-likelihood estimation (MLE) problems** in any scientific or engineering domain. The package does not assume a specific model type and can be applied to statistical models, mechanistic models, models defined through simulations, optimization problems, or arbitrary likelihood functions.

Typical application areas include (but are not limited to):

- Systems biology & Quantitative Systems Pharmacology,
- Engineering and control,
- Scientific Machine Learning and any field requiring identifiability and uncertainty analysis of MLE parameters.

## What problems does this package solve?

Profile likelihood methods provide insight into practical identifiability - how precisely model parameters (or predictions derived from them) are determined by the available data.
`LikelihoodProfiler.jl` offers a unified interface for:

- **Parameter profile likelihoods.** Profile the likelihood function to explore how well parameters are constrained by the data.
- **Confidence intervals for parameters.** Estimate confidence intervals for parameter values based on likelihood threshold to quantify the level of certainty in the parameter estimates.
- **Functional profile likelihoods.**  Profile arbitrary functions of the parameters (e.g., predictions, reparameterizations, etc).

These capabilities apply to any setting where an MLE objective (e.g., negative log-likelihood) can be evaluated.

## Installation

In Julia terminal run the following command:

```julia
import Pkg; Pkg.add("LikelihoodProfiler")
```

## Citation

Borisov I., Metelkin E. An Algorithm for Practical Identifiability Analysis and Confidence Intervals Evaluation Based on Constrained Optimization. 2018. October. ICSB2018. https://doi.org/10.13140/RG.2.2.18935.06563
