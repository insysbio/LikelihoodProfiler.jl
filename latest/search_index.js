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
    "text": "Currently supported Julia versions are 0.7, 1.0.To install the package from REPLjulia> import Pkg   # if you are on Julia 0.7, 1.0\n\njulia> Pkg.add(PackageSpec(url=\"https://github.com/insysbio/LikelihoodProfiler.jl.git\"))"
},

{
    "location": "index.html#Quick-start-1",
    "page": "Home",
    "title": "Quick start",
    "category": "section",
    "text": "using LikelihoodProfiler\n\n# testing profile function\nf(x) = 5.0 + (x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2\n\n# Calculate parameters intervals for first parameter component, x[1]\nres_1 = get_interval(\n  [3., 2., 2.1], # starting point\n  1,             # parameter component\n  f,             # profile function\n  :LIN_EXTRAPOL; # method\n  loss_crit = 9. # critical level\n  )\n#\n\n# Plot parameter profile x[1]\nusing Plots\nplotly()\nplot(res_1)(Image: plot_lin)"
},

{
    "location": "index.html#Intro-1",
    "page": "Home",
    "title": "Intro",
    "category": "section",
    "text": "The reliability and predictability of a kinetic systems biology (SB) and systems pharmacology (SP) model depends on the calibration of model parameters. Taking into account the lacking of data and the experimental variability the value of any parameter determined unambiguously. This results in characterization of parameter by \"confidence intervals\" or even \"non-identifiable\" parameters when the confidence interval is open. The package includes algorithms to perform practical identifiability analysis and evaluation confidence intervals using Profile Likelihood [2] which can be applied to complex SB/SP models. Results of the identifiability analysis can be used to qualify and calibrate parameters or to reduce the model."
},

{
    "location": "index.html#Objective-1",
    "page": "Home",
    "title": "Objective",
    "category": "section",
    "text": "The package introduces several original algorithms taking into account the following points:This algorithm does not assume that the likelihood function is differentiable at any point. This allows using derivation free and global methods of optimization which do not require the calculation of gradients.\nThe calculation of likelihood function is the most computationally expensive operation within the others. It becomes critical for large dynamic model used nowadays in systems biology.\nThe algorithm should calculate the confidence endpoint with the selected tolerance and must be optimal regarding likelihood function calls. The intermediate (non-endpoint) profile points is not important.\nThe algorithm should be stable for calculation both finite and infinite intervals. They should stop immediately (with the corresponding status) if parameter is not identifiable."
},

{
    "location": "index.html#Methods-overview-1",
    "page": "Home",
    "title": "Methods overview",
    "category": "section",
    "text": "This algorithms can be applied to complex kinetic models where function differentiability is not guaranteed and each likelihood estimation is computationally expensive.  The package introduces original \"one-pass\" algorithm: Confidence Intervals evaluation by Constrained Optimization 6 developed by the authors of this package. :CICO_ONE_PASS utilizes the Inequality-based Constrained Optimization [3-4] for efficient determination of confidence intervals and detection of “non-identifiable” parameters.  The \"multi-pass\" methods use extrapolation/interpolation of likelihood points to the critical level: linear (:LIN_EXTRAPOL) and quadratic (:QUADR_EXTRAPOL) approaches. They are also effective for both identifiable and non-identifiable parameters."
},

{
    "location": "index.html#References-1",
    "page": "Home",
    "title": "References",
    "category": "section",
    "text": "Wikipedia Identifiability_analysis\nKreutz C., et al. Profile Likelihood in Systems Biology. FEBS Journal 280(11), 2564-2571, 2013\nSteven G. Johnson, The NLopt nonlinear-optimization package, link\nAndrew R. Conn, Nicholas I. M. Gould, and Philippe L. Toint, \"A globally convergent augmented Lagrangian algorithm for optimization with general constraints and simple bounds,\" SIAM J. Numer. Anal. vol. 28, no. 2, p. 545-572 (1991)\nJulia: A Fresh Approach to Numerical Computing. Jeff Bezanson, Alan Edelman, Stefan Karpinski and Viral B. Shah (2017) SIAM Review, 59: 65–98\nBorisov I., Metelkin E. An Algorithm for Practical Identifiability Analysis and Confidence Intervals Evaluation Based on Constrained Optimization. 2018. October. ICSB2018. https://doi.org/10.13140/RG.2.2.18935.06563"
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
    "text": "Three methods implemented in LikelihoodProfiler package: :CICO_ONE_PASS,  :LIN_EXTRAPOL,  :QUADR_EXTRAPOL.All the methods are wrapped by get_endpoint interface where the argument method should be initiated by one of the possible values. get_endpoint provides also the additional mechanism for better optimization by transforming the parameters to log/logit scales."
},

{
    "location": "methods.html#:CICO_ONE_PASS-1",
    "page": "Methods",
    "title": ":CICO_ONE_PASS",
    "category": "section",
    "text": "The method uses the one-pass calculation of confidence interval endpoint, i.e. one optimization is required for single endpoint. It utilizes the Inequality-based Constrained Optimization for efficient determination of confidence intervals and detection of “non-identifiable” parameters.The method internally calls NLopt algorithm to build an objective function with LN_AUGLAG algorithm."
},

{
    "location": "methods.html#:LIN_EXTRAPOL-1",
    "page": "Methods",
    "title": ":LIN_EXTRAPOL",
    "category": "section",
    "text": "The method uses multi-pass approach creating profile likelihood function and evaluating next step as linear extrapolation: y = a*x + b."
},

{
    "location": "methods.html#:QUADR_EXTRAPOL-1",
    "page": "Methods",
    "title": ":QUADR_EXTRAPOL",
    "category": "section",
    "text": "The method uses multi-pass approach creating profile likelihood function and evaluating next step as quadratic extrapolation: y = x^2*a + x*b + c."
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
    "text": "LikelihoodProfiler.get_interval function outputs estimated confidence interval along with other data as LikelihoodProfiler.ParamInterval structure.LikelihoodProfiler provides a @recipe for Plots.jl to visualize confidence interval estimation and plot parameter profile based on LikelihoodProfiler.ParamInterval.using LikelihoodProfiler\n\n# Likelihood function\nf(x) = 5.0 + (x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2\n\n# Calculate parameters intervals for x[1], x[2], x[3]\nres = [\n    get_interval(\n        [3., 2., 2.1],\n        i,\n        f,\n        :CICO_ONE_PASS;\n        loss_crit = 9.\n    ) for i in 1:3]\n\n# Plot parameter profile x[2]\nusing Plots\nplotly()\nplot(res[2])(Image: )To make a smooth plot compute more profile points using LikelihoodProfiler.update_profile_points! which internally uses PlotUtils.adapted_gridupdate_profile_points!(res[2])\n\nplot(res[2])(Image: )"
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
    "text": "function get_endpoint(\n    theta_init::Vector{Float64},\n    theta_num::Int,\n    loss_func::Function,\n    method::Symbol,\n    direction::Symbol = :right;\n\n    loss_crit::Float64 = 0.0,\n    scale::Vector{Symbol} = fill(:direct, length(theta_init)),\n    theta_bounds::Vector{Tuple{Float64,Float64}} = unscaling.(\n        fill((-Inf, Inf), length(theta_init)),\n        scale\n        ),\n    scan_bound::Float64 = unscaling(\n        (direction==:left) ? -9.0 : 9.0,\n        scale[theta_num]\n        ),\n    scan_tol::Float64 = 1e-3,\n    loss_tol::Float64 = 1e-3,\n    local_alg::Symbol = :LN_NELDERMEAD,\n    kwargs...\n    )\n\nCalculates right or left endpoint of CI for parameter component. It is a wripper of get_right_endpoint functions for selection of direction and using different transformations for faster optimization.\n\nReturn\n\nEndPoint object storing confidence endpoint and profile points found on fly.\n\nArguments\n\ntheta_init: starting values of parameter vector theta. The starting values is not necessary to be the optimum values for loss_func but it the value of loss_func must be lower than loss_crit.\ntheta_num: number n of vector component to compute confidence interval theta^n.\nloss_func: loss function Lambdaleft(thetaright) the profile of which is analyzed. Usually we use log-likelihood for profile analysis in form Lambda( theta ) = - 2 lnleft( L(theta) right).\nmethod: computational method to evaluate interval endpoint. Currently the following methods are implemented: :CICO_ONE_PASS, :LIN_EXTRAPOL, :QUADR_EXTRAPOL.\ndirection: :right or :left endpoint to estimate.\n\nKeyword arguments\n\nsee get_interval\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.get_interval-Tuple{Array{Float64,1},Int64,Function,Symbol}",
    "page": "API",
    "title": "LikelihoodProfiler.get_interval",
    "category": "method",
    "text": "get_interval(\n    theta_init::Vector{Float64},\n    theta_num::Int,\n    loss_func::Function,\n    method::Symbol;\n\n    loss_crit::Float64 = 0.0,\n    scale::Vector{Symbol} = fill(:direct, length(theta_init)),\n    theta_bounds::Vector{Tuple{Float64,Float64}} = unscaling.(\n        fill((-Inf, Inf), length(theta_init)),\n        scale\n        ),\n    scan_bounds::Tuple{Float64,Float64} = unscaling.(\n        (-9.0, 9.0),\n        scale[theta_num]\n        ),\n    scan_tol::Float64 = 1e-3,\n    loss_tol::Float64 = 1e-3,\n    local_alg::Symbol = :LN_NELDERMEAD,\n    kwargs...\n    )\n\nComputes confidence interval for single component theta_num of parameter vector and loss_func according to loss_crit level.\n\nReturn\n\nParamInterval structure storing all input data and estimated confidence interval.\n\nArguments\n\ntheta_init: starting values of parameter vector theta. The starting values is not necessary to be the optimum values for loss_func but it the value of loss_func must be lower than loss_crit.\ntheta_num: number n of vector component to compute confidence interval theta^n.\nloss_func: loss function Lambdaleft(thetaright) the profile of which is analyzed. Usually we use log-likelihood for profile analysis in form Lambda( theta ) = - 2 lnleft( L(theta) right).\nmethod: computational method to evaluate interval endpoint. Currently the following methods are implemented: :CICO_ONE_PASS, :LIN_EXTRAPOL, :QUADR_EXTRAPOL.\n\nKeyword arguments\n\nloss_crit: critical level of loss function. The endpoint of CI for selected parameter is the value at which profile likelihood meets the value of loss_crit.\nscale: vector of scale transformations for each component. Possible values: :direct, :log, :logit. This option can make optimization much more faster, especially for wide theta_bounds. The default value is :direct (no transformation) for all components.\ntheta_bounds: vector of bounds for each component in format (left_border, right_border). This bounds define the ranges for possible parameter values. The defaults are the non-limited values taking into account the scale, i.e. (0 Inf) for :log scale.\nscan_bounds: vector of scan bound for theta_num component. It must be within the theta_bounds for the scanned component. The defaults are (-9 9) for transformed values, i.e. (1e-9 1e9) for :log scale.\nscan_tol: Absolute tolerance of scanned component (stop criterion).\nloss_tol: Absolute tolerance of loss_func at loss_crit (stop criterion). Restriction. Currently is not effective for :CICO_ONE_PASS methods because of limitation in LN_AUGLAG interface.\nlocal_alg: algorithm of optimization. Currently the local derivation free algorithms form NLOPT pack were tested. The methods: :LN_NELDERMEAD, :LN_COBYLA, :LN_PRAXIS show good results. Methods: :LN_BOBYQA, :LN_SBPLX, :LN_NEWUOA is not recommended.\nkwargs...: the additional keyword arguments passed to get_right_endpoint for specific method.\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.get_right_endpoint",
    "page": "API",
    "title": "LikelihoodProfiler.get_right_endpoint",
    "category": "function",
    "text": "function get_right_endpoint(\n    theta_init::Vector{Float64},\n    theta_num::Int,\n    loss_func::Function,\n    method::Val{:CICO_ONE_PASS};\n\n    theta_bounds::Vector{Tuple{Float64,Float64}} = fill(\n        (-Inf, Inf), length(theta_init)\n        ),\n    scan_bound::Float64 = 9.0,\n    scan_tol::Float64 = 1e-3,\n    loss_tol::Float64 = 0.,\n    local_alg::Symbol = :LN_NELDERMEAD,\n    kwargs...\n    )\n\nInterface for current and future methods for endpoint estimation.\n\nReturn\n\nTuple of three values:\n\nRight end point value: ::Float64.\nProfile points estimated on fly: ::Array{ ProfilePoint, 1}, see ProfilePoint.\nStatus of sulution: ::Symbol. One of values: :BORDER_FOUND_BY_SCAN_TOL, :SCAN_BOUND_REACHED.\n\nArguments\n\ntheta_init: starting values of parameter vector theta. The starting values is not necessary to be the optimum values for loss_func but it the value of loss_func must be lower than loss_crit.\ntheta_num: number n of vector component to compute confidence interval theta^n.\nloss_func: loss function the profile of which is analyzed, see get_interval. In this function loss crit is always equal 0 for code simplification.\nmethod: this value is always fixed. Implemented methods are: Val{:CICO_ONE_PASS}. It is implemented for easy switching between different implemented and future methods.\n\nKeyword arguments\n\ntheta_bound: vector of bounds for each component in format (left_bound, right_bound). This bounds define the ranges for possible parameter values.\nscan_bound: right scan bound for theta_num component. It must be within the theta_bounds for the scanned component.\nscan_tol: Absolute tolerance of scanned component (stop criterion).\nloss_tol: Absolute tolerance of loss_func at loss_crit (stop criterion). Restriction. Currently is not effective for :CICO_ONE_PASS methods because of limitation in LN_AUGLAG interface.\nlocal_alg: algorithm of optimization. Currently the local derivation free algorithms form NLOPT pack were tested. The methods: :LN_NELDERMEAD, :LN_COBYLA, :LN_PRAXIS show good results. Methods: :LN_BOBYQA, :LN_SBPLX, :LN_NEWUOA is not recommended.\nkwargs...: the additional keyword arguments passed to get_right_endpoint for specific method.\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.profile-Tuple{Array{Float64,1},Int64,Function}",
    "page": "API",
    "title": "LikelihoodProfiler.profile",
    "category": "method",
    "text": "function profile(\n    theta_init::Vector{Float64},\n    theta_num::Int,\n    loss_func::Function;\n\n    skip_optim::Bool = false,\n    theta_bounds::Vector{Tuple{Float64,Float64}} = fill((-Inf, Inf), length(theta_init)),\n    local_alg::Symbol = :LN_NELDERMEAD,\n    ftol_abs::Float64 = 1e-3,\n    maxeval::Int = 10^5,\n    kwargs... # currently not used\n    )\n\nIt generates the profile function based on loss_func. Used internally in methods :LIN_EXTRAPOL, :QUADR_EXTRAPOL.\n\nReturn\n\nReturns profile function for selected parameter component. Each call of the function starts optimization.\n\nArguments\n\ntheta_init: starting values of parameter vector theta.\ntheta_num: number n of vector component to create the profile.\nloss_func: loss function Lambdaleft(thetaright) the profile of which is analyzed. Usually we use log-likelihood for profile analysis in form Lambda( theta ) = - 2 lnleft( L(theta) right).\n\nKeyword arguments\n\nskip_optim : set true if you need marginal profile, i.e. profile without optimization. Default is false.\ntheta_bounds : vector of bounds for each component in format (left_border, right_border). This bounds define the ranges for possible parameter values.\nlocal_alg : algorithm of optimization. Currently the local derivation free algorithms form NLOPT pack were tested. The methods: :LN_NELDERMEAD, :LN_COBYLA, :LN_PRAXIS show good results. Methods: :LN_BOBYQA, :LN_SBPLX, :LN_NEWUOA is not recommended.\nftol_abs : absolute tolerance criterion for profile function.\nmaxeval : maximal number of loss_func calls to estimate profile point.\n\n\n\n\n\n"
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
    "text": "update_profile_points!(\n    pi::ParamInterval;\n    max_recursions::Int = 2\n    )\n\nRefines profile points to make your plot more smooth. Internally uses adapted_grid to compute additional profile points. See PlotUtils.adapted_grid.\n\nArguments\n\npi : ParamInterval structure to update.\nmax_recursions : how many times each interval is allowed to\n\nbe refined (default: 2).\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.EndPoint",
    "page": "API",
    "title": "LikelihoodProfiler.EndPoint",
    "category": "type",
    "text": "struct EndPoint\n    value::Float64\n    profilePoints::Array{ProfilePoint, 1}\n    status::Symbol\n    direction::Symbol\n    counter::Int\nend\n\nStructure storing end point for confidence interval.\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.ParamInterval",
    "page": "API",
    "title": "LikelihoodProfiler.ParamInterval",
    "category": "type",
    "text": "struct ParamInterval\n    input::ParamIntervalInput\n    loss_init::Float64\n    method::Symbol\n    result::Tuple{EndPoint, EndPoint}\nend\n\nStructure storing result of parameter interval calculation\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.ParamIntervalInput",
    "page": "API",
    "title": "LikelihoodProfiler.ParamIntervalInput",
    "category": "type",
    "text": "struct ParamIntervalInput\n    theta_init::Vector{Float64} # initial parameters vector\n    theta_num::Int # number of the parameter for analysis\n    loss_func::Function # loss function\n    loss_crit::Float64 # loss function maximum value, \"identifiability level\"\n    scale::Vector{Symbol}\n    theta_bounds::Vector{Tuple{Float64, Float64}} # search bounds for id parameter\n    scan_bounds::Tuple{Float64,Float64}\n    scan_tol::Float64 # fitting tolerance for local optimizer (default - 1e-3)\n    loss_tol::Float64 # constraints tolerance\n    local_alg::Symbol # local fitting algorithm (default - :LN_NELDERMEAD)\n    fitter_options::Any\nend\n\nStructure storing input data for parameter interval calculation\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.ProfilePoint",
    "page": "API",
    "title": "LikelihoodProfiler.ProfilePoint",
    "category": "type",
    "text": "struct ProfilePoint\n    value::Float64\n    loss::Float64\n    params::Array{Float64, 1}\n    ret::Symbol\n    counter::Union{Int, Nothing}\nend\n\nStructure storing one point from profile function.\n\n\n\n\n\n"
},

{
    "location": "api.html#LikelihoodProfiler.LikelihoodProfiler",
    "page": "API",
    "title": "LikelihoodProfiler.LikelihoodProfiler",
    "category": "module",
    "text": "Main module for LikelihoodProfiler.jl.\n\nTwo functions are exported from this module for public use:\n\nget_endpoint. Calculates endpoint of confidence interval.\nget_interval. Calculates confidence interval.\nprofile. Generates the profile function based on loss_func\nupdate_profile_points!. Updates interval by profile points.\n\n\n\n\n\n"
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
