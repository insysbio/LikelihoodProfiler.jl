# Methods

 Three methods implemented in **LikelihoodProfiler** package: `:CICO_ONE_PASS`,
 `:LIN_EXTRAPOL`,
 `:QUADR_EXTRAPOL`.

The main function for CI endpoints estimation is `get_interval` (see API section). Estimation method is defined by the keyword argument `method`. It supports one of the above values. Default is `:CICO_ONE_PASS`.

## :CICO\_ONE\_PASS

The method uses the one-pass calculation of confidence interval endpoint and  utilizes **Inequality-based Constrained Optimization**
for efficient determination of confidence intervals and detection of “non-identifiable” parameters.

 The method internally calls [`NLopt`](https://nlopt.readthedocs.io/en/latest/) algorithm to build an **augmented objective function** with [`LN_AUGLAG`](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#augmented-lagrangian-algorithm) algorithm. Optimization algorithm choice is described in Optimizers section.

## :LIN_EXTRAPOL

The method examines profile likelihood function by making steps in both directions from the optima. 
Each step is derived as a linear extrapolation: `y = a*x + b`.

## :QUADR_EXTRAPOL

The method examines profile likelihood function by making steps in both directions from the optima. 
Each step is derived as a quadratic extrapolation: `y = x^2*a + x*b + c`.
