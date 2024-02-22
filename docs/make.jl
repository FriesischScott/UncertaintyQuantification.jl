using DataFrames
using Documenter
using DocumenterCitations
using UncertaintyQuantification

DocMeta.setdocmeta!(
    UncertaintyQuantification,
    :DocTestSetup,
    :(using UncertaintyQuantification, DataFrames, DisplayAs, Random; Random.seed!(8128)),
)

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"))

makedocs(
    bib;
    modules=[UncertaintyQuantification],
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
    sitename="UncertaintyQuantification.jl",
    authors="Jasper Behrensdorf and Ander Gray",
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Getting Started" => "manual/gettingstarted.md",
            "Reliability Analysis" => "manual/reliability.md",
            "Metamodelling" => "manual/metamodels.md",
            "High Performance Computing" => "manual/hpc.md",
        ],
        "Examples" => [
            "Metamodels" => "examples/metamodels.md",
            "High Performance Computing" => "examples/hpc.md",
        ],
        "API" => [
            "Inputs" => "api/inputs.md",
            "Parameter" => "api/parameter.md",
            "RandomVariable" => "api/randomvariable.md",
            "ResponseSurface" => "api/responsesurface.md",
            "PolyharmonicSpline" => "api/polyharmonicspline.md",
            "Simulations" => "api/simulations.md",
            "SlurmInterface" => "api/slurm.md",
        ],
        "References" => "references.md",
    ],
    strict=:doctest,
)

deploydocs(;
    repo="github.com/FriesischScott/UncertaintyQuantification.jl", push_preview=true
)
