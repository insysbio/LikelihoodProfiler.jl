push!(LOAD_PATH,"../src/")

using Documenter, ParametersIdentification

makedocs(
    format = :html,
    sitename = "ParametersIdentification",
    pages = [
        "index.md",
        "basics.md"
    ]
)

deploydocs(
    repo   = "http://gitlab.insilicobio.ru/development/ParametersIdentification.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
