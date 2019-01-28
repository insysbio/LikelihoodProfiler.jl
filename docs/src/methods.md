# Methods

 Three methods implemented in **LikelihoodProfiler** package: `:CICO_ONE_PASS`,
 `:LIN_EXTRAPOL`,
 `:QUADR_EXTRAPOL`.

All the methods are wrapped by `get_endpoint` interface where the argument `method` should be initiated by one of the possible values. `get_endpoint` provides also the additional mechanism for better optimization by transforming the parameters to log/logit scales.

## :CICO\_ONE\_PASS

The method uses the one-pass calculation of confidence interval endpoint, i.e. one
optimization is required for single endpoint. It utilizes the **Inequality-based Constrained Optimization**
for efficient determination of confidence intervals and detection of “non-identifiable” parameters.

 The method internally calls [`NLopt`](https://nlopt.readthedocs.io/en/latest/) algorithm to build an **objective function** with [`LN_AUGLAG`](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#augmented-lagrangian-algorithm) algorithm.

## :LIN_EXTRAPOL

The method uses multi-pass approach creating profile likelihood function and evaluating
next step as linear extrapolation: `y = a*x + b`.

## :QUADR_EXTRAPOL

The method uses multi-pass approach creating profile likelihood function and evaluating
next step as quadratic extrapolation: `y = x^2*a + x*b + c`.
