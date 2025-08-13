---
title: "LikelihoodProfiler.jl: A Unified Julia Package for Practical Identifiability Analysis and Confidence Interval Estimation"
tags:
  - profile likelihood
  - identifiability analysis
  - parameter estimation
  - systems biology
  - quantitative systems pharmacology
  - mathematical modeling
authors:
  - name: "Ivan Borisov"
    orcid: 0000-0002-2887-1903
    affiliation: 1
  - name: "Evgeny Metelkin"
    orcid: 0000-0001-6612-5373
    affiliation: 1
  - name: "Aleksander Demin"
    affiliation: 1
affiliations:
  - name: "InSysBio LLC"
    index: 1
date: 12 August 2025
bibliography: paper.bib
---

## Summary

Practical identifiability addresses the critical question of how well a mechanistic model is determined by the available experimental data. In a typical model's calibration workflow practical identifiability implies prior structural identifiability, which studies and resolves the uncertainty within the model structure independently from the available data. However, in many cases, profile likelihood-based methods are also used as a proxy for structural identifiability analysis, particularly when the complexity of a model makes structural methods inapplicable or computationally prohibitive [@Wieland2021]. Moreover, profile likelihood techniques can be extended beyond parameter analysis to assess the identifiability of model states and predictions. This versatility makes profile likelihood analysis an essential component in the development and validation of models in Systems Biology and Quantitative Systems Pharmacology (QSP).

LikelihoodProfiler.jl is an open-source Julia package designed to perform profile likelihood-based identifiability analyses by offering a unified and extensible interface.

## Statement of Need

Despite the widespread use of profile likelihood methods in practical identifiability analysis, existing tools often lack a common interface and extensive ecosystem integration, limiting accessibility and reproducibility. Also different profile likelihood based methods are implemented in different software tools, which requires the user to switch between different languages and software interfaces. LikelihoodProfiler.jl addresses these limitations by providing:
- Unified interface to access multiple profiling methods.
- Compatibility with common modeling standards (Heta [@Metelkin2021], PEtab [@Persson2025], SBML).
- Integration with Julia’s SciML [@Rackauckas2017] for efficient computation and extensibility.

## Features and Methodologies

LikelihoodProfiler.jl supports the following methods for profiling likelihood functions and estimating confidence intervals:
1. **OptimizationProfiler** follows the classical approach, employing stepwise re-optimization of the likelihood function under parameter constraints. It is intuitive but may be computationally intensive.
2. **IntegrationProfiler** computes likelihood profiles by solving differential equations derived from optimization problems, yielding smooth profiles. This method is advantageous but demands accurate Hessian computations or approximations, posing challenges for large models.
3. **CICOProfiler** estimates confidence interval endpoints directly via constrained optimization, avoiding full-profile reconstruction, thus providing efficient confidence interval estimation for practical scenarios [@Borisov2020].

All methods leverage a CommonSolve interface [@Rackauckas2017] (`CommonSolve.solve()`), supporting global settings for parallelization and verbosity. Profiling results are directly visualizable using the Plots.jl package and exportable as DataFrames for further analysis.

## Demonstrative Example: JAK/STAT Signaling Pathway Model

LikelihoodProfiler.jl’s functionality and interfaces are demonstrated using the JAK/STAT signaling pathway model [@Boehm2014], which consists of 8 states and 9 parameters. The model and experimental data were sourced from the Benchmark-Models-PEtab repository and imported through the PEtab.jl interface.

```julia
using PEtab, Plots
petab_model = PEtabModel("Boehm_JProteomeRes2014.yaml")
petab_problem = PEtabODEProblem(petab_model)
```

To define the profile likelihood problem, we construct an OptimizationProblem instance and provide the optimal parameter values:

```julia
using Optimization, LikelihoodProfiler
optprob = OptimizationProblem(petab_problem)
plprob = PLProblem(optprob, get_x(petab_problem))
```

LikelihoodProfiler.jl offers a suite of methods for profiling likelihood functions and assessing practical identifiability. Each method includes several configurable options, such as optimizer or integrator selection, tolerances, and step size. The most straightforward method is OptimizationProfiler, which follows the classical approach of stepwise re-optimization of the likelihood function under a constraint on the parameter of interest.

```julia
alg1 = OptimizationProfiler(optimizer = Optimization.LBFGS(), stepper = LineSearchStep(; initial_step=0.07))
```

A more advanced method is the IntegrationProfiler, which computes likelihood profiles by solving a system of differential equations derived from the underlying optimization problem. This method requires a differential equation solver (integrator) to be specified.

```julia
using OrdinaryDiffEq
alg2 = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.07,), matrix_type = :hessian)
```

An alternative approach, implemented in CICOProfiler, estimates the confidence intervals (CI) endpoints directly—without reconstructing the full profile—by solving a constrained optimization problem [@Borisov2020]. This method is often more efficient when only the CI is required.

```julia
alg3 = CICOProfiler(optimizer = :LN_NELDERMEAD, scan_tol = 1e-10)
```

All profiling methods share a common `solve()` interface.

```julia
sol = solve(plprob, alg1)
```

All three profiling approaches yielded comparable confidence intervals, emphasizing the reliability and flexibility of LikelihoodProfiler.jl for different modeling scenarios.

Below are the profile likelihoods for the first three parameters of the JAK/STAT model, computed using the three methods:

*figure: Profile likelihoods for the first three parameters of the JAK/STAT model using OptimizationProfiler, IntegrationProfiler, and CICOProfiler*

All three methods reported similar CI for the JAK/STAT model, which can be accessed using the `get_endpoint()` function.

The optimal profiling method and settings depend on the complexity of the model and the goal of the analysis:
- OptimizationProfiler benefits from the choice of optimization algorithm (e.g., gradient-based or derivative-free) but may be computationally intensive.
- IntegrationProfiler provides smooth profile trajectories but requires Hessian computation or approximation, which may be challenging for large-scale models.
- CICOProfiler is often more efficient for CI estimation when the full profile is not needed.

## Implementation and Extensibility

All profiling methods benefit from the unified interface provided by LikelihoodProfiler.jl:
- Integration with SciML packages gives users access to a wide range of optimizers, differential equation solvers, and AD backends, enabling efficient profiling configurations.
- Compatibility with Heta, PEtab and SBML formats broadens the accessibility of the package across different modeling frameworks.
- `solve()` interface provided by `CommonSolve.jl` provides unified access various profiling methods 
- A common parallelization setup, controlled via the `parallel_type` argument in the `solve()` function, is supported across all methods and can significantly accelerate computations.
- The interface facilitates integration of new profiling methods and stepping algorithms.

Future work will include adding new methods of parameters, functions and predictions profiling and enabling adaptive switching between strategies.

## Availability

LikelihoodProfiler.jl is open-source and available at: https://github.com/insysbio/LikelihoodProfiler.jl

The package is registered in Julia and can be installed from the Julia REPL using:
```julia
import Pkg; Pkg.add("LikelihoodProfiler")
```

Tutorials and documentation are available at: https://insysbio.github.io/LikelihoodProfiler.jl/stable/

## Acknowledgments

## References