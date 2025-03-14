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
        "API" => "api.md",
    ],
)


deploydocs(
    repo   = "github.com/insysbio/LikelihoodProfiler.jl.git"
)
