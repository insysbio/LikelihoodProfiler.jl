
using Documenter, LikelihoodProfiler

makedocs(
    modules = [LikelihoodProfiler],
    sitename = "LikelihoodProfiler.jl",
    pages = [
    "Home" => "index.md",
    "Methods" => "methods.md",
    "Visualization" => "visualization.md",
    "API" => "api.md",
    ],
)


deploydocs(
    repo   = "github.com/insysbio/LikelihoodProfiler.jl.git"
)
