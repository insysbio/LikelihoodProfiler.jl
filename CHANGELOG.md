# Change log

# 1.3.0

- Drop `Optimization.jl` and use `OptimizationBase.jl`
- Examples and tests modified in order to support new `Optimization.jl` libs
- `reoptimize_init` kwarg fixed
- Docs and README updated to reflect JOSS
- SIR model computational issues fixed
- Gaussian model added to the Tutorials

# 1.2.1

- Fixing bug in profiling bounds idxs

# 1.2.0

### Breaking changes
- `idxs` moved to problem level
- `profile_range` replaced with `profile_lower/profile_upper`
### Other changes
- target constructors added
- functions profiling added

## 1.1.3

- Zenodo release
- joss script fixed

## 1.1.2

- JOSS paper added

## 1.1.1

- `CICOBase` version updated to `v0.6.0` to allow equal `scan_bounds` and `theta_bounds`

## 1.1.0

- `PLProblem` renamed to `ProfileLikelihoodProblem` and `profile` to `solve`
- `CommonSolve.solve` interface added
- `LineSearchStep` drafted

## 1.0.2

- parallel modes added
- docs updated

## 1.0.1

- Basic section added to the docs
- Support optprob with parameters
- minor fixes (issues)

## 1.0.0

- new LikelihoodProfiler interface
- `CICOProfiler` method moved to Ext
- `OptimizationProfiler` added
- `IntegrationProfiler` added