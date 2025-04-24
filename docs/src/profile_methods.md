
## Profile likelihood methods

LikelihoodProfiler provides a range of methods to profile likelihood functions and explore practical identifiability. The method should be provided as the second argument to the [`profile`](@ref) function.

### Optimization-based profiles

The method computes profiles for each parameter by iteratively changing the value of the parameter and re-optimizing the likelihood function with respect to all other parameters. 

```@docs; canonical=false
OptimizationProfiler
```

### Integration-based profiles

The method computes profiles for each parameter by integrating the differential equations system. 

```@docs; canonical=false
IntegrationProfiler
```

References:
1. Chen, J.-S. & Jennrich, R. I. Simple Accurate Approximation of Likelihood Profiles. Journal of Computational and Graphical Statistics 11, 714–732 (2002).
2. Chen, J.-S. & Jennrich, R. I. The Signed Root Deviance Profile and Confidence Intervals in Maximum Likelihood Analysis. Journal of the American Statistical Association 91, 993–998 (1996).

### Confidence Intervals by Constrained Optimization (CICO)

The method computes intersections (endpoints of the confidence interval (CI)) of the profile with the predefined confidence level (`threshold`) without restoring the exact trajectory of the profile. Requires using [CICOBase](https://github.com/insysbio/CICOBase.jl) package.

```@docs; canonical=false
CICOProfiler 
```

References:
1. Borisov, I. & Metelkin, E. Confidence intervals by constrained optimization—An algorithm and software package for practical identifiability analysis in systems biology. PLoS Comput Biol 16, e1008495 (2020).
2. Venzon, D. J. & Moolgavkar, S. H. A Method for Computing Profile-Likelihood-Based Confidence Intervals. Applied Statistics 37, 87 (1988).

