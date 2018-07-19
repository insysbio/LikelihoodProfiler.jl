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
    "text": "ParametersIdentification is a Julia package for identifiability analysis and confidence intervals evaluation."
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
    "location": "index.html#References:-1",
    "page": "Overview",
    "title": "References:",
    "category": "section",
    "text": "1.Kreutz C., et al. Profile Likelihood in Systems Biology. FEBS Journal 280(11), 2564-2571, 20132.Steven G. Johnson, The NLopt nonlinear-optimization package, http://ab-initio.mit.edu/nlopt3.Andrew R. Conn, Nicholas I. M. Gould, and Philippe L. Toint, \"A globally convergent augmented Lagrangian algorithm for optimization with general constraints and simple bounds,\" SIAM J. Numer. Anal. vol. 28, no. 2, p. 545-572 (1991)4.Julia: A Fresh Approach to Numerical Computing. Jeff Bezanson, Alan Edelman, Stefan Karpinski and Viral B. Shah (2017) SIAM Review, 59: 65–98"
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
    "location": "basics.html#ParametersIdentification.params_intervals",
    "page": "Basics",
    "title": "ParametersIdentification.params_intervals",
    "category": "function",
    "text": "params_intervals(init_params::Vector{Float64},\n        id::Int64,\n        maxf::Float64,\n        loss_func::Function,\n        logscale::Vector{Bool};\n        fit_alg::Symbol=:LN_AUGLAG,\n        local_alg::Symbol=:LN_SBPLX,\n        init_bounds::Vector{Float64} = [1e-9,1e9],\n        tol_val::Float64 = 1e-4,\n        solver::Symbol=:NLOPT)\n\nInput:\n\n    init_params - initial parameters vector\n    id - id of the parameter for analysis\n    maxf - loss function maximum value, \"identifiability level\"\n    loss_func - loss function\n    logscale - bool vector length(init_params) where true - log scale / false - direct scale\n    fit_alg - fitting algorithm (default - :LN_AUGLAG)\n    local_alg - local fitting algorithm (default - :LN_SBPLX)\n    init_bounds - search bounds (default - [1e-9,1e9])\n    tol_val - fitting tolerance (default - 1e-4)\n    solver - fitting solver (default - :NLOPT)\n\nReturn:\n\n    confidence intervals evaluation:\n    (interval, termination reason, numer of evaluations)\n\n\n\n"
},

{
    "location": "basics.html#ParametersIdentification.params_plot",
    "page": "Basics",
    "title": "ParametersIdentification.params_plot",
    "category": "function",
    "text": "params_plot(init_params::Vector{Float64},\n            id::Int64,\n            interval::Tuple{Float64,Float64},\n            loss_func::Function,\n            maxf::Float64;\n            fit_alg::Symbol=:LN_NELDERMEAD,\n            tol_val::Float64=1e-4)\n\nInput:\n\n    init_params - initial parameters vector\n    id - id of the parameter for analysis\n    interval - interval for plot\n    loss_func - loss function\n    maxf - loss function maximum value, \"identifiability level\"\n    fit_alg - fitting algorithm (default - :LN_NELDERMEAD)\n    tol_val - fitting tolerance (default - 1e-4)\n\nReturn:\n\n    parameter profile plot\n\n\n\n"
},

{
    "location": "basics.html#Functions-1",
    "page": "Basics",
    "title": "Functions",
    "category": "section",
    "text": "ParametersIdentification.params_intervalsParametersIdentification.params_plot"
},

]}
