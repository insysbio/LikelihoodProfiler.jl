
using Documenter, LikelihoodProfiler

makedocs(
    format = :html,
    build   = "",
    sitename = "LikelihoodProfiler",
    pages = [
        "index.md",
        "basics.md"
    ],
    repo = "https://github.com/insysbio/LikelihoodProfiler.jl/blob/master/docs/src/index.md"
)

deploydocs(
    repo   = "",
    julia = "0.6",
    osname = "windows",
    target = "build",
    deps   = nothing,
    make = nothing
)
