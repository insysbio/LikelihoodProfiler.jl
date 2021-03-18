
using LikelihoodProfiler
using Documenter

makedocs(
    modules = [LikelihoodProfiler],
    sitename = "LikelihoodProfiler Documentation",
    pages = [
        "Home" => "index.md",
        "Methods" => "methods.md",
        "Visualization" => "visualization.md",
        "Optimizers" => "optimizers.md",
        "API" => "api.md",
    ],
)


deploydocs(
    repo   = "github.com/insysbio/LikelihoodProfiler.jl.git"
)
