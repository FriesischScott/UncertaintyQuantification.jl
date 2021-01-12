using Documenter, UncertaintyQuantification

makedocs(
    modules = [UncertaintyQuantification],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename="UncertaintyQuantification.jl",
    authors = "Jasper Behrensdorf and Ander Gray",
    doctest=false,
    pages = [
        "Home" => "index.md",
        "Uncertainty Propagation" => "UncProp.md",
        "Sensitivity Analysis" => "Sensitivity.md",
        "Reliability Analysis" => "Reliability.md",
        "Examples" => "examples.md",
        "API" => "api.md"
        ]
    )



deploydocs(
    repo   = "github.com/FriesischScott/UncertaintyQuantification.jl",
)