# Mean of a Gaussian Distribution

This tutorial demonstrates the basic workflow of `LikelihoodProfiler.jl` using one of the simplest statistical models: estimating the mean of a Gaussian distribution with known variance.

## Model and Data

We assume we observe i.i.d. data:
$$
X_1, \ldots, X_n \sim \mathcal{N}(\mu, 1)
$$
and want to construct a confidence interval for the mean parameter μ. 

```@example gaussian-1
using LikelihoodProfiler, OptimizationLBFGSB, Distributions, Random
using Plots

Random.seed!(73612768)

n_obs = 10
data = rand(Normal(0, 1), n_obs)
```

## Objective Function (Negative Log-Likelihood)

With known variance $\sigma^2 = 1$, the log-likelihood for μ is:

$$
\ell(\mu) = \sum_{i=1}^{n} \log \mathcal{N}(x_i \mid \mu, 1)
$$
We minimize negative log-likelihood:
```@example gaussian-1
obj(x, p) = -sum(logpdf.(Normal(x[1], 1.0), data))
```

## Maximum Likelihood Estimate

```@example gaussian-1
x0 = [1.0]   # initial guess

optf = OptimizationFunction(obj, AutoForwardDiff())
optprob = OptimizationProblem(optf, x0)
sol = solve(optprob, LBFGSB())

optpars = sol.u   # MLE of μ
```

For a Gaussian with known variance, this should be very close to the sample mean.

## Profile Likelihood for the Mean

We now build the `ProfileLikelihoodProblem`, define profiling method and solve the problem.
In this model there is one parameter, so the profile is one-dimensional.
```@example gaussian-1
plprob = ProfileLikelihoodProblem( optprob, optpars; 
    profile_lower = -3.0, profile_upper = 3.0)

method = OptimizationProfiler(optimizer = LBFGSB(),
    stepper = FixedStep(; initial_step = 0.01))

sol = solve(plprob, method)
plot(sol, size=(800,300), margins=5Plots.mm)
```

## Interpreting the plot

- The curve shows how the negative log-likelihood increases as μ moves away from its MLE.
- The horizontal line marks the confidence threshold for the chosen level (default 95%).
- The intersection points give the profile likelihood confidence interval.

## Extracting Confidence Interval Endpoints

`ProfileLikelihoodSolution` stores the computed profile curves together with confidence interval endpoints and identification retcodes, which indicate whether each parameter is practically identifiable. These values can be accessed directly:
 ```@example gaussian-1
retcodes(sol)
endpoints(sol)
```