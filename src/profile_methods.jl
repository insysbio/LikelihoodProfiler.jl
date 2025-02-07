abstract type AbstractProfilerMethod end

############################## OptimizationProfiler ##############################

"""
    OptimizationProfiler{S, opType, optsType}

A profiler method that uses stepwise re-optimization to profile the likelihood function.

### Fields

- `stepper::S`: The algorithm used to compute the next profile point.
- `optimizer::opType`: The optimizer used for the optimization process.
- `optimizer_opts::optsType`: Options for the optimizer. Defaults to `NamedTuple()`.

### Example 

```julia
using Optimization
profiler = OptimizationProfiler(; optimizer = Optimization.LBFGS(), optimizer_opts = (reltol=1e-4,), stepper = FixedStep())
```
"""
Base.@kwdef struct OptimizationProfiler{S, opType, optsType} <: AbstractProfilerMethod
  stepper::S
  optimizer::opType
  optimizer_opts::optsType = NamedTuple()
end
#!SciMLBase.allowsbounds(optalg) && error("Optimization alg used with `OptimizationProfiler` must support bounds.")

get_stepper(op::OptimizationProfiler) = op.stepper
get_optimizer(op::OptimizationProfiler) = op.optimizer
get_optimizer_opts(op::OptimizationProfiler) = op.optimizer_opts

############################## IntegrationProfiler ##############################
"""
    IntegrationProfiler{opType, optsType, DEAlg, DEOpts}

A profiler method that uses integration of differential equations system to profile the likelihood function.

### Fields

- `reoptimize::Bool`: Indicates whether to re-optimize during profiling. Defaults to `false`.
- `optimizer::opType`: The optimizer used for the optimization process. Defaults to `nothing`.
- `optimizer_opts::optsType`: Options for the optimizer. Defaults to `nothing`.
- `integrator::DEAlg`: The differential equation algorithm used for integration.
- `integrator_opts::DEOpts`: Options for the differential equation solver.
- `matrix_type::Symbol`: The type of matrix to be used for the Hessian approximation. Defaults to `:hessian`.
- `gamma::Float64`: Correction factor. Defaults to `1.0`.

### Example

```julia
using OrdinaryDiffEq
profiler = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.3,), matrix_type = :hessian)
```
"""
Base.@kwdef struct IntegrationProfiler{opType, optsType, DEAlg, DEOpts} <: AbstractProfilerMethod
  reoptimize::Bool = false
  optimizer::opType = nothing
  optimizer_opts::optsType = NamedTuple()
  integrator::DEAlg
  integrator_opts::DEOpts = NamedTuple()
  matrix_type::Symbol = :hessian
  gamma::Float64 = 1.0
end

get_integrator(ip::IntegrationProfiler) = ip.integrator
get_integrator_opts(ip::IntegrationProfiler) = ip.integrator_opts
get_matrix_type(ip::IntegrationProfiler) = ip.matrix_type
get_gamma(ip::IntegrationProfiler) = ip.gamma

############################## CICOProfiler ##############################

"""
    CICOProfiler

Confidence Intervals by Constrained Optimization (CICO) method to find the intersections of the likelihood function with the threshold.

### Fields

- `optimizer::Symbol`: The optimizer used for the optimization process. Defaults to NLopt `:LN_NELDERMEAD`.
- `scan_tol::Float64`: The tolerance for the endpoints scan. Defaults to `1e-3`.

### Example

```julia
profiler = CICOProfiler(optimizer = :LN_NELDERMEAD, scan_tol = 1e-3)
```
"""
Base.@kwdef struct CICOProfiler <: AbstractProfilerMethod
  optimizer::Symbol = :LN_NELDERMEAD
  scan_tol::Float64 = 1e-3
end

get_optimizer(cp::CICOProfiler) = cp.optimizer
get_scan_tol(cp::CICOProfiler) = cp.scan_tol