var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "LikelihoodProfiler is a Julia package for identifiability analysis and confidence intervals evaluation."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Julia download page.Currently supported Julia versions are 0.7, 1.0To install the package from REPLjulia> import Pkg   # if you are on Julia 0.7, 1.0\n\njulia> Pkg.add(PackageSpec(url=\"https://github.com/insysbio/LikelihoodProfiler.jl.git\"))\n\njulia> using LikelihoodProfiler"
},

{
    "location": "index.html#Quick-start-1",
    "page": "Home",
    "title": "Quick start",
    "category": "section",
    "text": "using LikelihoodProfiler\n\n# Likelihood function\nf(x) = 5.0 + (x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2\n\n# Calculate parameters intervals for x[1], x[2], x[3]\nres = [\n    get_interval(\n        [3., 2., 2.1],\n        i,\n        f,\n        :LIN_EXTRAPOL;\n        loss_crit = 9.\n    ) for i in 1:3]\n\n# Plot parameter profile x[1]\nusing Plots\nplotly()\nplot(res[1])(Image: )"
},

{
    "location": "index.html#Objective-1",
    "page": "Home",
    "title": "Objective",
    "category": "section",
    "text": "The reliability and predictability of a kinetic systems biology (SB) model depends on the calibration of model parameters. Experimental data can be insufficient to determine all the parameters unambiguously. This results in “non-identifiable” parameters and parameters identifiable within confidence intervals. The package includes algorithms for parameters identification using Profile Likelihood [1] method which can be applied to complex SB models. Results of the identifiability analysis can be used to qualify and calibrate parameters or to reduce the model."
},

{
    "location": "index.html#Methods-Overview-1",
    "page": "Home",
    "title": "Methods Overview",
    "category": "section",
    "text": "This packages provides a number of algorithms for identifiability analysis and confidence intervals evaluation by Profile Likelihood method. Along with linear extrapolation (:LIN_EXTRAPOL) and quadratic extrapolation (:QUADR_EXTRAPOL) the package introduces Confidence Intervals evaluation by Constrained Optimization method (:CICO_ONE_PASS) developed by the authors of this package. :CICO_ONE_PASS utilizes the Inequality-based Constrained Optimization [2, 3] for efficient determination of confidence intervals and detection of “non-identifiable” parameters. This algorithm does not assume that the likelihood function is differentiable or can be calculated for any given parameters set. This algorithm can be applied to complex kinetic models where function differentiability is not guaranteed and each likelihood estimation is computationally expensive.   The package includes tools for parameters identifiability analysis, confidence intervals evaluation and results visualization."
},

{
    "location": "index.html#References-1",
    "page": "Home",
    "title": "References",
    "category": "section",
    "text": "Kreutz C., et al. Profile Likelihood in Systems Biology. FEBS Journal 280(11), 2564-2571, 2013\nSteven G. Johnson, The NLopt nonlinear-optimization package, link\nAndrew R. Conn, Nicholas I. M. Gould, and Philippe L. Toint, \"A globally convergent augmented Lagrangian algorithm for optimization with general constraints and simple bounds,\" SIAM J. Numer. Anal. vol. 28, no. 2, p. 545-572 (1991)\nJulia: A Fresh Approach to Numerical Computing. Jeff Bezanson, Alan Edelman, Stefan Karpinski and Viral B. Shah (2017) SIAM Review, 59: 65–98"
},

{
    "location": "methods.html#",
    "page": "Methods",
    "title": "Methods",
    "category": "page",
    "text": ""
},

{
    "location": "methods.html#Methods-1",
    "page": "Methods",
    "title": "Methods",
    "category": "section",
    "text": "Three methods implemented in LikelihoodProfiler package: :CICO_ONE_PASS,  :LIN_EXTRAPOL, :QUADR_EXTRAPOL."
},

{
    "location": "methods.html#:CICO_ONE_PASS-1",
    "page": "Methods",
    "title": ":CICO_ONE_PASS",
    "category": "section",
    "text": "The method uses the one-pass calculation of confidence interval endpoint, i.e. one optimization is required for single endpoint. It utilizes the Inequality-based Constrained Optimization for efficient determination of confidence intervals and detection of “non-identifiable” parameters. The method internally calls NLopt algorithm to build an objective function with LN_AUGLAG algorithm. "
},

{
    "location": "methods.html#:LIN_EXTRAPOL-1",
    "page": "Methods",
    "title": ":LIN_EXTRAPOL",
    "category": "section",
    "text": "The method uses multi-pass approach creating profile likelihood function and evaluating next step as linear extrapolation: y = ax + b."
},

{
    "location": "methods.html#:QUADR_EXTRAPOL-1",
    "page": "Methods",
    "title": ":QUADR_EXTRAPOL",
    "category": "section",
    "text": "The method uses multi-pass approach creating profile likelihood function and evaluating next step as quadratic extrapolation: y = x^2a + xb + c."
},

{
    "location": "methods.html#Methods-comparison-1",
    "page": "Methods",
    "title": "Methods comparison",
    "category": "section",
    "text": "The next results are generated automatically based on current version."
},

{
    "location": "visualization.html#",
    "page": "Visualization",
    "title": "Visualization",
    "category": "page",
    "text": ""
},

{
    "location": "visualization.html#Visualization-1",
    "page": "Visualization",
    "title": "Visualization",
    "category": "section",
    "text": "LikelihoodProfiler.get_interval function outputs estimated confidence interval along with other data as LikelihoodProfiler.ParamInterval structure.LikelihoodProfiler provides a @recipe for Plots.jl to visualize confidence interval estimation and plot parameter profile based on LikelihoodProfiler.ParamInterval.using LikelihoodProfiler\n\n# Likelihood function\nf(x) = 5.0 + (x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2\n\n# Calculate parameters intervals for x[1], x[2], x[3]\nres = [\n    get_interval(\n        [3., 2., 2.1],\n        i,\n        f,\n        :CICO_ONE_PASS;\n        loss_crit = 9.\n    ) for i in 1:3]\n\n# Plot parameter profile x[1]\nusing Plots\nplotly()\nplot(res[2])(Image: )To make a smooth plot compute more profile points using LikelihoodProfiler.update_profile_points! which internally uses PlotUtils.adapted_gridupdate_profile_points!(res[2])\n\nplot(res[2])(Image: )"
},

{
    "location": "api.html#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api.html#LikelihoodProfiler.get_endpoint",
    "page": "API",
    "title": "LikelihoodProfiler.get_endpoint",
    "category": "function",
    "text": "function get_endpoint(\n    theta_init::Vector{Float64},\n    theta_num::Int,\n    loss_func::Function,\n    method::Symbol,\n    direction::Symbol = :right;\n\n    loss_crit::Float64 = 0.0,\n    scale::Vector{Symbol} = fill(:direct, length(theta_init)),\n    theta_bounds::Vector{Tuple{Float64,Float64}} = unscaling.(\n        fill((-Inf, Inf), length(theta_init)),\n        scale\n        ),\n    scan_bound::Float64 = unscaling(\n        (direction==:left) ? -9.0 : 9.0,\n        scale[theta_num]\n        ),\n    scan_tol::Float64 = 1e-3,\n    loss_tol::Float64 = 1e-3,\n    local_alg::Symbol = :LN_NELDERMEAD,\n    kwargs...\n    )\n\nCalculates right or left endpoint of CI for parameter component. It is a wripper of get_right_endpoint functions for selection of direction and using different transformations for faster optimization.\n\nReturn\n\nEndPoint object storing confidence endpoint and profile points found on fly.\n\nArguments\n\ntheta_init: starting values of parameter vector theta. The starting values is not necessary to be the optimum values for loss_func but it the value of loss_func must be lower than loss_crit.\ntheta_num: number n of vector component to compute confidence interval theta^n.\nloss_func: loss function Lambdaleft(thetaright) the profile of which is analyzed. Usually we use log-likelihood for profile analysis in form Lambda( theta ) = - 2 lnleft( L(theta) right)\nmethod: computational method to evaluate interval endpoint. Currently the following methods are implemented: :CICO_ONE_PASS.\ndirection: :right or :left endpoint to estimate.\n\nKeyword arguments\n\nsee get_interval\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.get_interval-Tuple{Array{Float64,1},Int64,Function,Symbol}",
    "page": "API",
    "title": "LikelihoodProfiler.get_interval",
    "category": "method",
    "text": "get_interval(\n    theta_init::Vector{Float64},\n    theta_num::Int,\n    loss_func::Function,\n    method::Symbol;\n\n    loss_crit::Float64 = 0.0,\n    scale::Vector{Symbol} = fill(:direct, length(theta_init)),\n    theta_bounds::Vector{Tuple{Float64,Float64}} = unscaling.(\n        fill((-Inf, Inf), length(theta_init)),\n        scale\n        ),\n    scan_bounds::Tuple{Float64,Float64} = unscaling.(\n        (-9.0, 9.0),\n        scale[theta_num]\n        ),\n    scan_tol::Float64 = 1e-3,\n    loss_tol::Float64 = 1e-3,\n    local_alg::Symbol = :LN_NELDERMEAD,\n    kwargs...\n    )\n\nComputes confidence interval for single component theta_num of parameter vector and loss_func according to loss_crit level.\n\nReturn\n\nParamInterval structure storing all input data and estimated confidence interval.\n\nArguments\n\ntheta_init: starting values of parameter vector theta. The starting values is not necessary to be the optimum values for loss_func but it the value of loss_func must be lower than loss_crit.\ntheta_num: number n of vector component to compute confidence interval theta^n.\nloss_func: loss function Lambdaleft(thetaright) the profile of which is analyzed. Usually we use log-likelihood for profile analysis in form Lambda( theta ) = - 2 lnleft( L(theta) right)\nmethod: computational method to evaluate interval endpoint. Currently the following methods are implemented: :CICO_ONE_PASS.\n\nKeyword arguments\n\nloss_crit: critical level of loss function.\nscale: vector of scale transformations for each component. Possible values: :direct, :log, :logit. This option can make optimization much more faster, especially for wide theta_bounds.\ntheta_bounds: vector of bounds for each component in format (left_border, right_border). The values outside the bound will be ignored.\nscan_bounds: vector of scan bound for theta_num component. It must be within the theta_bounds for the scanned component.\nscan_tol: Abolute tolerance of scanned component (stop criterion).\nloss_tol: experimental. Required tolerance of loss_func at loss_crit.\nlocal_alg: algorithm of optimization. Currently the local derivation free algorithms form NLOPT pack were tested. The methods: :LN_NELDERMEAD, :LN_COBYLA, :LN_PRAXIS show good results. Methods: :LN_BOBYQA, :LN_SBPLX, :LN_NEWUOA is not recommended.\nkwargs...: experimental. other keyword arguments passed to optimization method.\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.get_right_endpoint",
    "page": "API",
    "title": "LikelihoodProfiler.get_right_endpoint",
    "category": "function",
    "text": "function get_right_endpoint(\n    theta_init::Vector{Float64},\n    theta_num::Int,\n    loss_func::Function,\n    method::Val{:CICO_ONE_PASS};\n\n    theta_bounds::Vector{Tuple{Float64,Float64}} = fill(\n        (-Inf, Inf), length(theta_init)\n        ),\n    scan_bound::Float64 = 9.0,\n    scan_tol::Float64 = 1e-3,\n    loss_tol::Float64 = 0.,\n    local_alg::Symbol = :LN_NELDERMEAD,\n    kwargs...\n    )\n\nInterface for current and future methods for endpoint estimation.\n\nReturn\n\nTuple of three values:\n\nRight end point value: ::Float64.\nProfile points estimated on fly: ::Array{ ProfilePoint, 1}, see ProfilePoint.\nStatus of sulution: ::Symbol. One of values: :BORDER_FOUND_BY_SCAN_TOL, :SCAN_BOUND_REACHED.\n\nArguments\n\ntheta_init: starting values of parameter vector theta. The starting values is not necessary to be the optimum values for loss_func but it the value of loss_func must be lower than loss_crit.\ntheta_num: number n of vector component to compute confidence interval theta^n.\nloss_func: loss function the profile of which is analyzed, see get_interval. In this function loss crit is always equal 0 for code simplification.\nmethod: this value is always fixed. Implemented methods are: Val{:CICO_ONE_PASS}. It is implemented for easy switching between different implemented and future methods.\n\nKeyword arguments\n\ntheta_bound: vector of bounds for each component in format (left_bound, right_bound). The values outside the bound will be ignored.\nscan_bound: right scan bound for theta_num component. It must be within the theta_bounds for the scanned component.\nscan_tol: Abolute tolerance of scanned component (stop criterion).\nloss_tol: experimental. Required tolerance of loss_func at loss_crit.\nlocal_alg: algorithm of optimization. Currently the local derivation free algorithms form NLOPT pack were tested. The methods: :LN_NELDERMEAD, :LN_COBYLA, :LN_PRAXIS show good results. Methods: :LN_BOBYQA, :LN_SBPLX, :LN_NEWUOA is not recommended.\nkwargs...: experimental. other keyword arguments passed to optimization method. max_iter: maximal number of loss_func calls. ftol_abs: absolute tolerance of parameter vector components.\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.profile-Tuple{Array{Float64,1},Int64,Function}",
    "page": "API",
    "title": "LikelihoodProfiler.profile",
    "category": "method",
    "text": "function profile(\n    theta_init::Vector{Float64},\n    theta_num::Int,\n    loss_func::Function;\n\n    skip_optim::Bool = false,\n    theta_bounds::Vector{Tuple{Float64,Float64}} = fill((-Inf, Inf), length(theta_init)),\n    local_alg::Symbol = :LN_NELDERMEAD,\n    ftol_abs::Float64 = 1e-3,\n    maxeval::Int = 10^5,\n    kwargs...\n    )\n\nReturns profile function for selected parameter component.\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.scaling",
    "page": "API",
    "title": "LikelihoodProfiler.scaling",
    "category": "function",
    "text": "scaling(x::Float64, scale::Symbol = :direct)\n\nTransforms values from specific scale to range [-Inf, Inf] based on option.\n\nReturn\n\nTransformed value.\n\nArguments\n\nx: input value.\nscale: transformation type: :direct, :log, :logit.\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.unscaling",
    "page": "API",
    "title": "LikelihoodProfiler.unscaling",
    "category": "function",
    "text": "unscaling(x::Float64, scale::Symbol = :direct)\n\nTransforms values from [-Inf, Inf] to specific scale based on option. Inverse function for scaling.\n\nReturn\n\nTransformed value.\n\nArguments\n\nx: input value.\nscale: transformation type: :direct, :log, :logit.\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.update_profile_points!-Tuple{ParamInterval}",
    "page": "API",
    "title": "LikelihoodProfiler.update_profile_points!",
    "category": "method",
    "text": "update_profile_points!(pi::ParamInterval)\n\nRefines profile points to make your plot more smooth. Internally uses adapted_grid to compute additional profile points. See PlotUtils.adapted_grid\n\nArguments\n\nmax_recursions: how many times each interval is allowed to\n\nbe refined (default: 2).\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.EndPoint",
    "page": "API",
    "title": "LikelihoodProfiler.EndPoint",
    "category": "type",
    "text": "struct EndPoint\n    value::Float64\n    profilePoints::Array{ProfilePoint, 1}\n    status::Symbol\n    direction::Symbol\n    counter::Int\nend\n\nEnd point storage.\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.ParamInterval",
    "page": "API",
    "title": "LikelihoodProfiler.ParamInterval",
    "category": "type",
    "text": "Structure storing result of parameter interval calculation\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.ParamIntervalInput",
    "page": "API",
    "title": "LikelihoodProfiler.ParamIntervalInput",
    "category": "type",
    "text": "Structure storing input data for parameter interval calculation\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.ProfilePoint",
    "page": "API",
    "title": "LikelihoodProfiler.ProfilePoint",
    "category": "type",
    "text": "Structure storing one point from profile function\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.LikelihoodProfiler",
    "page": "API",
    "title": "LikelihoodProfiler.LikelihoodProfiler",
    "category": "module",
    "text": "Main module for LikelihoodProfiler.jl.\n\nTwo functions are exported from this module for public use:\n\nget_endpoint. Calculates endpoint of confidence interval.\nget_interval. Calculates confidence interval.\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.adapted_grid2-Tuple{Any,Tuple{Real,Real}}",
    "page": "API",
    "title": "LikelihoodProfiler.adapted_grid2",
    "category": "method",
    "text": "adapted_grid2(f, minmax::Tuple{Number, Number}; max_recursions = 7)\n\nComputes a grid x on the interval [minmax[1], minmax[2]] so that plot(f, x) gives a smooth \"nice\" plot. The method used is to create an initial uniform grid (21 points) and refine intervals where the second derivative is approximated to be large. When an interval becomes \"straight enough\" it is no longer divided. Functions are never evaluated exactly at the end points of the intervals. The parameter max_recusions computes how many times each interval is allowed to be refined.\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.logistic10-Tuple{Float64}",
    "page": "API",
    "title": "LikelihoodProfiler.logistic10",
    "category": "method",
    "text": "logistic10(x::Float64)\n\nFunction transforming interval [-Inf, Inf] to [0,1] using logistic transformation. Inverse function for logit10.\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.logit10-Tuple{Float64}",
    "page": "API",
    "title": "LikelihoodProfiler.logit10",
    "category": "method",
    "text": "logit10(x::Float64)\n\nFunction transforming interval [0,1] to [-Inf, Inf] using logit transformation.\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.plot2-Tuple{Array{ProfilePoint,1}}",
    "page": "API",
    "title": "LikelihoodProfiler.plot2",
    "category": "method",
    "text": "Experimental. Plot array of ProfilePoint.\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.plot2-Tuple{EndPoint}",
    "page": "API",
    "title": "LikelihoodProfiler.plot2",
    "category": "method",
    "text": "Experimental. Plot EndPoint.\n\n\n\n\n\n"
},

{
    "location": "api.html#RecipesBase.apply_recipe-Tuple{Dict{Symbol,Any},ParamInterval}",
    "page": "API",
    "title": "RecipesBase.apply_recipe",
    "category": "method",
    "text": "using Plots\nplotly()\nplot(pi::ParamInterval)\n\nPlots profile L(theta) for parameter theta_num, identifiability level, identifiability interval. Use update_profile_points!(pi::ProfileInterval) function to refine profile points and make your plot more smooth\n\n\n\n\n\n"
},

{
    "location": "api.html#API-references-1",
    "page": "API",
    "title": "API references",
    "category": "section",
    "text": "The package exports the following functions for parameters identifiability analysis, confidence intervals evaluation and results visualization.Modules = [LikelihoodProfiler]\nOrder   = [:function, :type, :module]"
},

]}
