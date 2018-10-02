push!(LOAD_PATH,"../src/")

using Documenter, LikelihoodProfiler

makedocs(
    format = :html,
    sitename = "LikelihoodProfiler",
    pages = [
        "index.md",
        "basics.md"
    ]
)

deploydocs(
    repo   = "gitlab.insilicobio.ru/development/LikelihoodProfiler.git",
    julia = "0.6",
    osname = "windows",
    target = "build",
    deps   = not
)
