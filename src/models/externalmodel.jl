struct ExternalModel <: UQModel
    sourcedir::String
    sources::Vector{String}
    extractors::Vector{Extractor}
    solver::Solver
    workdir::String
    extras::Vector{String}
    formats::Dict{Symbol,String}
    cleanup::Bool

    function ExternalModel(
        sourcedir::String,
        sources::Union{String,Vector{String}},
        extractors::Union{Extractor,Vector{Extractor}},
        solver::Solver,
        workdir::String,
        extras::Union{String,Vector{String}},
        formats::Dict{Symbol,String},
        cleanup::Bool,
    )
        sources, extractors, extras = wrap.([sources, extractors, extras])
        return new(
            sourcedir, sources, extractors, solver, workdir, extras, formats, cleanup
        )
    end
end

function ExternalModel(
    sourcedir::String,
    sources::Union{String,Vector{String}},
    extractors::Union{Extractor,Vector{Extractor}},
    solver::Solver;
    workdir::String=tempname(),
    extras::Union{String,Vector{String}}=String[],
    formats::Dict{Symbol,String}=Dict{Symbol,String}(),
    cleanup::Bool=false,
)
    return ExternalModel(
        sourcedir, sources, extractors, solver, workdir, extras, formats, cleanup
    )
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
            if isempty(file)
                continue
            end
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

function formatinputs(row::DataFrameRow, formats::Dict{Symbol,String})
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
