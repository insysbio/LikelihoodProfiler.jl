using LikelihoodProfiler
using Documenter

makedocs(
    modules = [LikelihoodProfiler],
    sitename = "LikelihoodProfiler Documentation",
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Extended Tutorials" => [
            "Taxol model" => "case_studies/taxol.md",
            "Bachmann" => "case_studies/bachmann.md",
            "Large model 1" => "case_studies/large_model1.md",
            "Large model 2" => "case_studies/large_model2.md",
        ],
        "API" => "api.md",
    ],
)


deploydocs(
    repo   = "github.com/insysbio/LikelihoodProfiler.jl.git"
)
