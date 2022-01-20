using Documenter, UncertaintyQuantification

DocMeta.setdocmeta!(
    UncertaintyQuantification,
    :DocTestSetup,
    :(using UncertaintyQuantification, Random; Random.seed!(8128));
    recursive=true,
)

makedocs(;
    modules=[UncertaintyQuantification],
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
    sitename="UncertaintyQuantification.jl",
    authors="Jasper Behrensdorf and Ander Gray",
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Uncertainty Propagation" => "manual/uncertainty.md",
            "Sensitivity Analysis" => "manual/sensitivity.md",
            "Reliability Analysis" => "manual/reliability.md",
        ],
        "API" => [
            "Inputs" => "api/inputs.md",
            "Parameter" => "api/parameter.md",
            "RandomVariable" => "api/randomvariable.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/FriesischScott/UncertaintyQuantification.jl", push_preview=true
)
