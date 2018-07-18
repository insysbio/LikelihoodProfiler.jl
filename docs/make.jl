push!(LOAD_PATH,"../src/")

using Documenter, parameters_identification

makedocs(
    format = :html,
    sitename = "ParametersIdentification",
    pages = [
        "index.md"
    ]
)

deploydocs(
    repo   = "http://gitlab.insilicobio.ru/development/ParametersIdentification.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
