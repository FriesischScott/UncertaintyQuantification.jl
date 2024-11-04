using DataFrames
using Documenter
using DocumenterCitations
using UncertaintyQuantification

format = if !isempty(ARGS) && ARGS[1] == "vite"
    using DocumenterVitepress

    MarkdownVitepress(;
        repo="https://github.com/FriesischScott/UncertaintyQuantification.jl",
        deploy_url="https://friesischscott.github.io/UncertaintyQuantification.jl",
    )
else
    Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true")
end

DocMeta.setdocmeta!(
    UncertaintyQuantification,
    :DocTestSetup,
    :(using UncertaintyQuantification, DataFrames, DisplayAs, Random; Random.seed!(8128)),
)

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"))

makedocs(;
    modules=[UncertaintyQuantification],
    plugins=[bib],
    format=format,
    sitename="UncertaintyQuantification.jl",
    authors="Jasper Behrensdorf and Ander Gray",
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Introduction" => "manual/introduction.md",
            "Getting Started" => "manual/gettingstarted.md",
            "Reliability Analysis" => "manual/reliability.md",
            "Metamodelling" => "manual/metamodels.md",
            "Simulations" => "manual/simulations.md",
            "Bayesian Updating" => "manual/bayesianupdating.md",
            "Parallelisation" => "manual/parallelisation.md",
            "High Performance Computing" => "manual/hpc.md",
        ],
        "Examples" => [
            "External Models" => "examples/external.md",
            "Metamodels" => "examples/metamodels.md",
            "Bayesian Updating" => "examples/bayesianupdating.md",
            "High Performance Computing" => "examples/hpc.md",
        ],
        "Benchmarks" => ["Subset Simulation" => "benchmarks/subset.md"],
        "API" => [
            "Inputs" => "api/inputs.md",
            "Models" => "api/models.md",
            "ResponseSurface" => "api/responsesurface.md",
            "PolyharmonicSpline" => "api/polyharmonicspline.md",
            "Simulations" => "api/simulations.md",
            "Bayesian Updating" => "api/bayesianupdating.md",
            "SlurmInterface" => "api/slurm.md",
        ],
        "References" => "references.md",
    ],
    warnonly=false,
)

deploydocs(;
    repo="github.com/FriesischScott/UncertaintyQuantification.jl", push_preview=true
)
