

################################ SELECTORS ################################

get_optimizer(method::AbstractProfilerMethod) = nothing
get_integrator(method::AbstractProfilerMethod) = nothing
hasoptimizer(method::AbstractProfilerMethod) = !isnothing(get_optimizer(method))
hasintegrator(method::AbstractProfilerMethod) = !isnothing(get_integrator(method))

############################## OptimizationProfiler ##############################

"""
    OptimizationProfiler{S, opType, optsType}

A profiler method that uses stepwise re-optimization to profile the likelihood function.

### Fields

- `stepper::S`: The algorithm used to compute the next profile point. Supported steppers include:
    - `FixedStep`: Proposes a constant step size in the profiling direction (Default). 
    - `LineSearchStep`: Uses a line search to adaptively determine the step size in the direction which is chosen by the `direction` keyword argument.
- `optimizer::opType`: The optimizer used for the optimization process.
- `optimizer_opts::optsType`: Options for the optimizer. Defaults to `NamedTuple()`.

### Stepping Options

The `stepper` argument controls how the next profile point is chosen. For example:

- `stepper = FixedStep(initial_step=0.1)`: Use a constant step size of 0.1.
- `stepper = LineSearchStep(direction=:Secant, linesearch=InterpolationLineSearch())`: Use a line search with secant direction.

See the documentation for each stepper type (e.g., `?FixedStep`, `?LineSearchStep`) for more details and customization options.

### Example

```julia
using Optimization
profiler = OptimizationProfiler(; optimizer = Optimization.LBFGS(), optimizer_opts = (reltol=1e-4,))
```
"""
Base.@kwdef struct OptimizationProfiler{S, opType, optsType} <: AbstractProfilerMethod
  stepper::S = FixedStep()
  optimizer::opType
  optimizer_opts::optsType = NamedTuple()
end

get_optimizer(op::OptimizationProfiler) = op.optimizer
get_optimizer_opts(op::OptimizationProfiler) = op.optimizer_opts
get_stepper(op::OptimizationProfiler) = op.stepper

############################## IntegrationProfiler ##############################
"""
    IntegrationProfiler{opType, optsType, DEAlg, DEOpts}

A profiler method that uses integration of differential equations system to profile the likelihood function.

### Fields

- `reoptimize::Bool`: Indicates whether to re-optimization after each step of the `integrator`. Defaults to `false`.
- `optimizer::opType`: The optimizer used for the optimization process. Defaults to `nothing`.
- `optimizer_opts::optsType`: Options for the optimizer. Defaults to `NamedTuple()`.
- `integrator::DEAlg`: The differential equation algorithm used for integration.
- `integrator_opts::DEOpts`: Options for the differential equation solver. Defaults to `NamedTuple()`.
- `matrix_type::Symbol`: The type of matrix to be used for the Hessian approximation. Possible options are: `:hessian`, `:identity`. Defaults to `:hessian`.
- `gamma::Float64`: Correction factor used in integration if full hessian is not computed (e.g. `matrix_type = :identity`). Defaults to `1.0`.

### Example

```julia
using OrdinaryDiffEq
profiler = IntegrationProfiler(integrator = Tsit5(), integrator_opts = (dtmax=0.3,), matrix_type = :hessian)
```
"""
struct IntegrationProfiler{opType, optsType, DEAlg, DEOpts} <: AbstractProfilerMethod
  reoptimize::Bool
  optimizer::opType
  optimizer_opts::optsType
  integrator::DEAlg
  integrator_opts::DEOpts
  matrix_type::Symbol
  gamma::Float64
end

function IntegrationProfiler(; reoptimize::Bool=false,
                               optimizer=nothing,
                               optimizer_opts=NamedTuple(),
                               integrator,                     # required
                               integrator_opts=NamedTuple(),
                               matrix_type::Symbol=:hessian,
                               gamma::Real=1.0)

  # matrix_type
  (matrix_type === :identity || matrix_type === :fisher || matrix_type === :hessian) ||
      throw(ArgumentError("`matrix_type` must be one of :identity, :fisher, :hessian (got $matrix_type)."))

  # reoptimize requires an optimizer
  reoptimize && optimizer === nothing &&
      throw(ArgumentError("`reoptimize=true` requires an `optimizer` to be provided."))

  # gamma sanity
  gamma_val = float(gamma)
  gamma_val >= 0 || throw(ArgumentError("Correction factor `gamma` must be non-negative (got $gamma_val)."))

  return IntegrationProfiler{typeof(optimizer), typeof(optimizer_opts),
                             typeof(integrator), typeof(integrator_opts)}(
           reoptimize, optimizer, optimizer_opts, integrator, integrator_opts,
           matrix_type, gamma_val)
end

get_optimizer(ip::IntegrationProfiler) = ip.optimizer
get_optimizer_opts(ip::IntegrationProfiler) = ip.optimizer_opts
get_integrator(ip::IntegrationProfiler) = ip.integrator
get_integrator_opts(ip::IntegrationProfiler) = ip.integrator_opts
get_matrix_type(ip::IntegrationProfiler) = ip.matrix_type
get_gamma(ip::IntegrationProfiler) = ip.gamma
############################## CICOProfiler ##############################

"""
    CICOProfiler

Confidence Intervals by Constrained Optimization (CICO) method to find the intersections of the likelihood function with the threshold.
See [CICOBase docs](https://github.com/insysbio/CICOBase.jl) for more details. Requires `using CICOBase`.

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