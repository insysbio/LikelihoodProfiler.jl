## FIMProfiler proposal

`FIMProfiler` is intended as a fast, local (asymptotic) confidence-interval method based on the Fisher information matrix (FIM), complementary to full profile-likelihood methods.

### Scope

- `FIMProfiler` should return **Wald-type confidence intervals** around the optimum.
- It should be explicitly documented as an **approximation** (local quadratic log-likelihood assumption).
- It should not try to reconstruct profile trajectories.

### Proposed user API

```julia
profiler = FIMProfiler(; 
    inversion=:cholesky,
    clamp_to_bounds=true,
)

sol = solve(plprob, profiler)
```

### Fields

- `resolve_fim(plprob, method::FIMProfiler)` is exported for users who want direct access to the matrix used by `FIMProfiler`.
- `inversion::Symbol`
  - `:cholesky`, `:pinv`, or `:svd`.
- `clamp_to_bounds::Bool`
  - clip intervals to profile bounds.

### Algorithm

Assume objective is negative log-likelihood and optimum is `θ̂`.

1. Compute/obtain Hessian `H = ∇² objective(θ̂)`.
2. Symmetrize: `H = (H + H') / 2`.
3. Estimate covariance: `Σ = inv(F)` (or pseudoinverse strategy).
4. Standard errors: `seᵢ = sqrt(max(Σᵢᵢ, 0))`.
5. Quantile proxy from problem threshold: `z = sqrt(plprob.threshold)`.
6. CI for profiled index `i`: `[θ̂ᵢ - z*seᵢ, θ̂ᵢ + z*seᵢ]`.
7. Optional: clamp to user bounds (`profile_lower/profile_upper`).

### Integration into LikelihoodProfiler internals

- Keep public interface unchanged: `solve(plprob, method; kwargs...)`.
- Add specialized internal dispatch:

```julia
__solve(plprob::ProfileLikelihoodProblem, method::FIMProfiler; kwargs...)
```

This specialization should bypass branch construction `[(idx, dir)]` and parallel profile stepping.

### Returned solution

To preserve downstream plotting/access patterns, construct synthetic 2-point profiles:

- per parameter `i`, create x-grid `[lbᵢ, ubᵢ]` from the estimated CI endpoints;
- assign local quadratic objective values relative to optimum,
  `Δℓ(x) = ((x - θ̂ᵢ)^2) / (2*seᵢ^2)`.

This keeps the return type compatible while clearly representing a local approximation.

### Numerical diagnostics (important)

`FIMProfiler` should attach diagnostics and warnings:

- condition number of `H`;
- smallest eigenvalue;
- whether pseudoinverse was used;
- which parameters have near-zero/negative variance estimates;
- if intervals were clipped by bounds.

### Behavior recommendations

- Ignore `parallel_type` (or warn once): no branch-wise parallelization needed.
- Throw a clear error when Hessian is unavailable and cannot be generated.
- For constrained/transformed parameters, compute CIs in optimization space first, then map back carefully (future extension).

### Why this design

- InformationGeometry uses Fisher-metric/Hessian-based geometry as a local metric and also offers non-adaptive profile modes that estimate profile domain from inverse Fisher information; this supports using FIM for fast local uncertainty estimates.
- PEtab workflows expose gradients/Hessians and are commonly used with profile-likelihood tooling, making observed-information CIs a practical fast pre-screen before full profiling.
- In OED workflows (e.g. Corleone family tooling), Fisher-information based uncertainty/optimality criteria are standard, so this profiler aligns with experiment-design practice.

### Planned phased implementation

1. **Phase 1 (MVP)**: observed-information CIs from Hessian at optimum.
2. **Phase 2**: robust inversion modes and richer diagnostics.
3. **Phase 3**: expected-FIM hooks for model-specific backends.
4. **Phase 4**: transformed-parameter and prediction CI support.

### Should FIM be stored in `ProfileLikelihoodProblem`?

Current design keeps `ProfileLikelihoodProblem` minimal: `FIMProfiler` uses the Hessian pathway from `OptimizationProblem` (thus reusing user backend choices). If PEtab integration requires custom FIM/FIM! hooks later, this can be added as an explicit extension point.
