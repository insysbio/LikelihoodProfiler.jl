using LikelihoodProfiler
using Documenter

makedocs(
    modules = [LikelihoodProfiler],
    sitename = "LikelihoodProfiler.jl",
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Extended Tutorials" => [
            "Taxol model" => "case_studies/taxol.md"
        ],
        "Basics" => [
            "Problem interface" => "problem_interface.md",
            "Profile likelihood methods" => "profile_methods.md",
            "Solution interface" => "solution_interface.md",
        ],
        "API" => "api.md",
    ],
)


deploydocs(
    repo   = "github.com/insysbio/LikelihoodProfiler.jl.git"
)
