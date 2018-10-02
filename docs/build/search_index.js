var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Overview-1",
    "page": "Overview",
    "title": "Overview",
    "category": "section",
    "text": "LikelihoodProfiler is a Julia package for identifiability analysis and confidence intervals evaluation."
},

{
    "location": "index.html#Installation-1",
    "page": "Overview",
    "title": "Installation",
    "category": "section",
    "text": "Julia download page. Currently supported Julia versions are 0.6, 0.7To install the package from REPLjulia> Pkg.add(PackageSpec(url=\"http://gitlab.insilicobio.ru/development/LLikelihoodProfiler.git\"))\r\n\r\njulia> using LikelihoodProfiler"
},

{
    "location": "index.html#Objective-1",
    "page": "Overview",
    "title": "Objective",
    "category": "section",
    "text": "The reliability and predictability of a kinetic systems biology (SB) model depends on the calibration of model parameters. Experimental data can be insufficient to determine all the parameters unambiguously. This results in “non-identifiable” parameters and parameters identifiable within confidence intervals. The proposed algorithm is a practical implementation of Profile Likelihood [1] method for parameters identification which can be applied to complex SB models. The results of this algorithm can be used to qualify and calibrate parameters or to reduce the model."
},

{
    "location": "index.html#Algorithm-1",
    "page": "Overview",
    "title": "Algorithm",
    "category": "section",
    "text": "The proposed algorithm for Profile Likelihood method addresses the disadvantages and restrictions of the root-finding algorithms with regard to the above problem and utilizes the Inequality-based Constrained Optimization [2, 3] for efficient determination of confidence intervals and detection of “non-identifiable” parameters. This algorithm does not assume that the likelihood function is differentiable or can be calculated for any given parameters set. This algorithm can be applied to complex kinetic models where function differentiability is not guaranteed and each likelihood estimation is computationally expensive.  The algorithm was tested for the set of kinetic models and it is distributed as a software package based on Julia Programming Language [4]. The package includes tools for parameters identifiability analysis, confidence intervals evaluation and results visualization."
},

{
    "location": "index.html#References-1",
    "page": "Overview",
    "title": "References",
    "category": "section",
    "text": "Kreutz C., et al. Profile Likelihood in Systems Biology. FEBS Journal 280(11), 2564-2571, 2013\nSteven G. Johnson, The NLopt nonlinear-optimization package, link\nAndrew R. Conn, Nicholas I. M. Gould, and Philippe L. Toint, \"A globally convergent augmented Lagrangian algorithm for optimization with general constraints and simple bounds,\" SIAM J. Numer. Anal. vol. 28, no. 2, p. 545-572 (1991)\nJulia: A Fresh Approach to Numerical Computing. Jeff Bezanson, Alan Edelman, Stefan Karpinski and Viral B. Shah (2017) SIAM Review, 59: 65–98"
},

{
    "location": "basics.html#",
    "page": "Basics",
    "title": "Basics",
    "category": "page",
    "text": ""
},

{
    "location": "basics.html#Basics-1",
    "page": "Basics",
    "title": "Basics",
    "category": "section",
    "text": "The package exports the following functions for parameters identifiability analysis, confidence intervals evaluation and results visualization."
},

{
    "location": "basics.html#LikelihoodProfiler.params_intervals",
    "page": "Basics",
    "title": "LikelihoodProfiler.params_intervals",
    "category": "function",
    "text": "params_intervals(init_params::Vector{Float64}, id::Int64,\nloss_crit::Float64, loss_func::Function; <keyword arguments>)\n\nComputes confidence interval for id parameter of init_params vector and loss_func according to loss_crit confidence level.\n\nReturns ParamInterval structure storing all input data and estimated confidence interval.\n\nArguments\n\nmethod::Symbol: computational method (:ONE_PASS,:D2D_PLE).\nlogscale_all::Bool: set logscale for all parameters to true / false (default false).\nlogscale::Vector{Bool}: set logscale for each parameter (default false for all parameters).\nscan_bound::Vector{Float64}: search bounds for id parameter (default [-9.,9.]).\nlocal_alg::Symbol: fitting algorithm (default :LN_NELDERMEAD).\nbounds::Vector{Vector{Float64}}: bound constraints for all parameters (default [-Inf,Inf]).\nmax_iter::Int64: maximum loss_func evaluations (default 10^5).\nptol::Float64: fitting tolerance for optimizer (default 1e-3).\nlosstol::Float64: constraints tolerance (default 1e-3).\n\n\n\n\n\n"
},

{
    "location": "basics.html#LikelihoodProfiler.params_plot",
    "page": "Basics",
    "title": "LikelihoodProfiler.params_plot",
    "category": "function",
    "text": "params_plot(params::Vector{Float64}, id::Int64, loss_func::Function,\ninterval::Tuple{Float64,Float64}; <keyword arguments>)\n\nComputes adapted_grid for loss_func and id parameter values from the interval. See also: Plots.adapted_grid\n\nArguments\n\nfit_alg::Symbol: fitting algorithm (default :LN_NELDERMEAD).\nbounds::Vector{Vector{Float64}}: bound constraints for all parameters (default [-Inf,Inf]).\ntol::Float64: fitting tolerance (default ftol_abs = 1e-3).\nmax_recursions::Int64: how many times each interval is allowed to be refined (default 2).\n\n\n\n\n\n"
},

{
    "location": "basics.html#Functions-1",
    "page": "Basics",
    "title": "Functions",
    "category": "section",
    "text": "LikelihoodProfiler.params_intervalsLikelihoodProfiler.params_plot"
},

]}
