abstract type AbstractStochasticProcess <: UQInput end

struct SpectralRepresentation <: AbstractStochasticProcess
    psd::AbstractPowerSpectralDensity
    time::AbstractVector{<:Real}
    name::Symbol
end

function sample(sr::SpectralRepresentation, n::Integer=1)
    return DataFrame(
        sr.name => map(collect, eachcol(rand(Uniform(0, 2π), (length(sr.psd.ω), n))))
    )
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

function to_standard_normal_space!(sr::SpectralRepresentation, df::DataFrame)
    df[!, sr.name] = map(p -> quantile.(Normal(), cdf.(Uniform(0, 2π), p)), df[:, sr.name])
    return nothing
end

function to_physical_space!(sr::SpectralRepresentation, df::DataFrame)
    df[!, sr.name] = map(p -> quantile.(Uniform(0, 2π), cdf.(Normal(), p)), df[:, sr.name])
    return nothing
end
