using LikelihoodProfiler
using Documenter

makedocs(
    modules = [LikelihoodProfiler],
    sitename = "LikelihoodProfiler.jl",
    pages = [
        "LikelihoodProfiler.jl: Unified Interface to Profile Likelihood Methods" => "index.md",
        "Getting Started" => "rosenbrock.md",
        "Tutorials" => [
            "Taxol model" => "case_studies/taxol.md"
        ],
        "Basics" => [
            "Problem interface" => "problem_interface.md",
            "Profile likelihood methods" => "profile_methods.md",
            "Solution interface" => "solution_interface.md",
            "Parallel execution" => "parallel_modes.md"
        ],
        "API" => "api.md",
    ],
)


deploydocs(
    repo   = "github.com/insysbio/LikelihoodProfiler.jl.git"
)
