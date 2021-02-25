struct ExternalModel <: UQModel
    sourcedir::String
    sources::Array{String,1}
    extras::Array{String,1}
    formats::Dict{Symbol,FormatSpec}
    workdir::String
    extractors::Array{Extractor,1}
    solver::Solver
end

function evaluate!(m::ExternalModel, df::DataFrame)

    datetime = Dates.format(now(), "YYYY-mm-dd-HH-MM-SS")

    n = size(df, 1)
    digits = ndigits(n)

    results = pmap(1:n) do i
        path = joinpath(m.workdir, datetime, "sample-$(lpad(i, digits, "0"))")
        mkpath(path)

        row = formatinputs(df[i, :], m.formats)

        for file ∈ m.sources
            tokens = Mustache.load(joinpath(m.sourcedir, file))

            open(joinpath(path, file), "w") do io
                render(io, tokens, row)
            end
        end

        for file ∈ m.extras
            cp(joinpath(m.sourcedir, file), joinpath(path, file))
        end

        run(m.solver, path)

        return map(e -> e.f(path), m.extractors)
    end

    results = hcat(results...) |> transpose
    vars = names(m.extractors)

    for (i, name) ∈ enumerate(names(m.extractors))
        df[!, name] = results[:, i]
    end
end

function formatinputs(row::DataFrameRow, formats::Dict{Symbol,FormatSpec})
    names = propertynames(row)
    values = []
    for symbol ∈ names
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