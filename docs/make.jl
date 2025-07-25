using DataFrames
using Documenter
using UncertaintyQuantification

if isempty(ARGS) || ARGS[1] != "pdf"
    using DocumenterCitations
    bib = CitationBibliography(joinpath(@__DIR__, "references.bib"))
end

format = if !isempty(ARGS) && ARGS[1] == "vite"
    using DocumenterVitepress
    MarkdownVitepress(;
        repo="github.com/FriesischScott/UncertaintyQuantification.jl",
        devbranch="master",
        devurl="dev",
    )
elseif !isempty(ARGS) && ARGS[1] == "pdf"
    include("./TypstWriter.jl")
    using .TypstWriter

    TypstWriter.Typst(;
        platform="native",
        bib=joinpath(@__DIR__, "references.bib"),
        version=pkgversion(UncertaintyQuantification),
    )
else
    Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true")
end

DocMeta.setdocmeta!(
    UncertaintyQuantification,
    :DocTestSetup,
    :(using UncertaintyQuantification, DataFrames, DisplayAs, Random; Random.seed!(8128)),
)

makedocs(;
    modules=[UncertaintyQuantification],
    plugins=isempty(ARGS) || ARGS[1] != "pdf" ? [bib] : [],
    format=format,
    sitename="UncertaintyQuantification.jl",
    authors="Jasper Behrensdorf and Ander Gray",
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Introduction" => "manual/introduction.md",
            "Getting Started" => "manual/gettingstarted.md",
            "Kernel Density Estimation" => "manual/kde.md",
            "Reliability Analysis" => "manual/reliability.md",
            "Metamodelling" => "manual/metamodels.md",
            "Simulations" => "manual/simulations.md",
            "Bayesian Updating" => "manual/bayesianupdating.md",
            "Parallelisation" => "manual/parallelisation.md",
            "Stochastic Dynamics" => "manual/dynamics.md",
            "High Performance Computing" => "manual/hpc.md",
        ],
        "Examples" => [
            "External Models" => "examples/external.md",
            "Metamodels" => "examples/metamodels.md",
            "Bayesian Updating" => "examples/bayesianupdating.md",
            "Stochastic Dynamics" => "examples/dynamics.md",
            "High Performance Computing" => "examples/hpc.md",
        ],
        "Benchmarks" => ["Subset Simulation" => "benchmarks/subset.md"],
        "API" => [
            "Inputs" => "api/inputs.md",
            "Models" => "api/models.md",
            "Reliability" => "api/reliability.md",
            "ResponseSurface" => "api/responsesurface.md",
            "PolyharmonicSpline" => "api/polyharmonicspline.md",
            "Simulations" => "api/simulations.md",
            "Bayesian Updating" => "api/bayesianupdating.md",
            "Power Spectral Density Functions" => "api/psd.md",
            "Stochastic Processes (Spectral Representation)" => "api/spectralrepresentation.md",
            "SlurmInterface" => "api/slurm.md",
        ],
        "References" => "references.md",
    ],
    warnonly=false,
    draft=false,
    source="src",
    build="build",
)

deploydocs(;
    repo="github.com/FriesischScott/UncertaintyQuantification.jl", push_preview=true
)
