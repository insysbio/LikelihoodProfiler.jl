# LikelihoodProfiler.jl: Unified Interface to Profile Likelihood Methods

`LikelihoodProfiler.jl` is a [Julia](https://julialang.org/downloads/) package for practical identifiability analysis and confidence intervals estimation using the profile likelihood approach. The package provides a unified interface for various profile likelihood methods, including optimization-based **OptimizationProfiler** and integration-based profiles **IntegrationProfiler**, CI endpoints search **CICOProfiler**, and more.

## Who is this package for?

`LikelihoodProfiler.jl` is intended for researchers and practitioners working with **maximum-likelihood estimation (MLE) problems** in any scientific or engineering domain. The package does not assume a specific model type and can be applied to statistical models, mechanistic models, models defined through simulations, optimization problems, or arbitrary likelihood functions. Applicable to a broad class of MLE problems, profile likelihood is especially valuable for complex nonlinear models — such as high-dimensional mechanistic models — where standard asymptotic confidence intervals can be unreliable.

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

## Related packages

Other implementations of the profile likelihood approach in Julia include:
- [ProfileLikelihood.jl](https://github.com/DanielVandH/ProfileLikelihood.jl) implements fixed-step optimization-based profiles and supports bivariate profile likelihoods.
- [InformationalGeometry.jl](https://github.com/RafaelArutjunjan/InformationGeometry.jl) implements various methods to study likelihood functions (including profile likelihood) using the tools of differential geometry.

There are also well-known profile likelihood implementations in other languages, namely: [Data2Dynamics](https://github.com/Data2Dynamics/d2d), [dMod](https://github.com/dkaschek/dMod/), [pyPESTO](https://github.com/ICB-DCM/pyPESTO), [sbioparametersci](https://www.mathworks.com/help/simbio/ref/sbioparameterci.html)

## Citation

Borisov I., Metelkin E. An Algorithm for Practical Identifiability Analysis and Confidence Intervals Evaluation Based on Constrained Optimization. 2018. October. ICSB2018. https://doi.org/10.13140/RG.2.2.18935.06563
