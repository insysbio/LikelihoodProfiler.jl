
using Documenter, LikelihoodProfiler

makedocs(
    format = :html,
    build   = "",
    sitename = "LikelihoodProfiler",
    pages = [
        "index.md",
        "basics.md"
    ]
)

deploydocs(
    repo   = "",
    julia = "0.6",
    osname = "windows",
    target = "build",
    deps   = not
)
