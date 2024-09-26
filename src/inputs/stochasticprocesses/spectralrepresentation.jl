abstract type AbstractStochasticProcess <: RandomUQInput end

struct SpectralRepresentation <: AbstractStochasticProcess
    psd::AbstractPowerSpectralDensity
    time::AbstractVector{<:Real}
    name::Symbol
end

function sample(sr::SpectralRepresentation, n::Integer=1)
    return DataFrame(names(sr) .=> eachcol(rand(Uniform(0, 2π), (n, length(sr.psd.ω)))))
end

function evaluate(sr::SpectralRepresentation, ϕ::AbstractVector{<:Real})
    S = evaluate(sr.psd)

    Δω = sr.psd.ω[2] - sr.psd.ω[1]
    A = sqrt.(2 * S * Δω)

    A[sr.psd.ω .== 0.0] .= 0.0

    return reduce(
        +,
        [sqrt(2) * an .* cos.(ωn .* sr.time .+ ϕn) for (an, ωn, ϕn) in zip(A, sr.psd.ω, ϕ)],
    )
end

function (sr::SpectralRepresentation)(ϕ::AbstractVector{<:Real})
    return evaluate(sr, ϕ)
end

function evaluate(sr::SpectralRepresentation, row::DataFrameRow)
    ϕ = collect(row[[Symbol("$(sr.name)_$i") for i in 1:length(sr.psd.ω)]])

    return evaluate(sr, ϕ)
end

function (sr::SpectralRepresentation)(row::DataFrameRow)
    return evaluate(sr, row)
end

function to_standard_normal_space!(sr::SpectralRepresentation, df::DataFrame)
    for v in names(sr)
        df[!, v] = quantile.(Normal(), cdf.(Uniform(0, 2π), df[:, v]))
    end
    return nothing
end

function to_physical_space!(sr::SpectralRepresentation, df::DataFrame)
    for v in names(sr)
        df[!, v] = quantile.(Uniform(0, 2π), cdf.(Normal(), df[:, v]))
    end
    return nothing
end

function dimensions(sr::SpectralRepresentation)
    return length(sr.psd.ω)
end

function names(sr::SpectralRepresentation)
    return [Symbol("$(sr.name)_$i") for i in 1:dimensions(sr)]
end
