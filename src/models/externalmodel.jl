struct ExternalModel <: UQModel
    sourcedir::String
    sources::Vector{String}
    extras::Vector{String}
    formats::Dict{Symbol,FormatSpec}
    workdir::String
    extractors::Vector{Extractor}
    solver::Solver
    cleanup::Bool

    function ExternalModel(
        sourcedir::String,
        sources::Vector{String},
        extras::Vector{String},
        formats::Dict{Symbol,FormatSpec},
        workdir::String,
        extractors::Vector{Extractor},
        solver::Solver,
        cleanup::Bool=false,
    )
        return new(
            sourcedir, sources, extras, formats, workdir, extractors, solver, cleanup
        )
    end
end

function evaluate!(
    m::ExternalModel,
    df::DataFrame;
    datetime::String=Dates.format(now(), "YYYY-mm-dd-HH-MM-SS"),
)
    n = size(df, 1)
    digits = ndigits(n)

    results = pmap(1:n) do i
        path = joinpath(m.workdir, datetime, "sample-$(lpad(i, digits, "0"))")
        mkpath(path)

        row = formatinputs(df[i, :], m.formats)

        for file in m.sources
            tokens = Mustache.load(joinpath(m.sourcedir, file))

            open(joinpath(path, file), "w") do io
                render(io, tokens, row)
            end
        end

        for file in m.extras
            cp(joinpath(m.sourcedir, file), joinpath(path, file))
        end

        run(m.solver, path)

        result = map(e -> e.f(path), m.extractors)
        if m.cleanup
            rm(path; recursive=true)
        end
        return result
    end

    results = hcat(results...)

    for (i, name) in enumerate(names(m.extractors))
        df[!, name] = results[i, :]
    end
end

function formatinputs(row::DataFrameRow, formats::Dict{Symbol,FormatSpec})
    names = propertynames(row)
    values = []
    for symbol in names
        if haskey(formats, symbol)
            push!(values, fmt(formats[symbol], row[symbol]))
        elseif haskey(formats, :*)
            push!(values, fmt(formats[:*], row[symbol]))
        else
            push!(values, row[symbol])
        end
    end
    return (; zip(names, values)...)
end
