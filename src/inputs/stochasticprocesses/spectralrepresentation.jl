"""
    AbstractStochasticProcess <: RandomUQInput

An abstract type representing a stochastic process in the context of uncertainty quantification. Concrete types inheriting from this abstract type should implement specific models for stochastic processes.

# Subtypes
- `SpectralRepresentationProcess`
- `OtherStochasticProcess`

# Description
Stochastic processes are mathematical models used to describe systems that evolve over time with inherent randomness. This abstract type serves as a base for defining various stochastic process models, which can be used in simulations and analyses within the uncertainty quantification framework.

"""

abstract type AbstractStochasticProcess <: RandomUQInput end

"""
    SpectralRepresentation(psd, time, name)

Alternative constructor

```julia
    SpectralRepresentation(psd, time, name)  # `islog` = true
```

"""

struct SpectralRepresentation <: AbstractStochasticProcess
    psd::AbstractPowerSpectralDensity
    time::AbstractVector{<:Real}
    ωt::AbstractMatrix{<:Real}
    Δω::Real
    A::AbstractVector{<:Real}
    name::Symbol
    ϕnames::Vector{Symbol}
end

function SpectralRepresentation(
    psd::AbstractPowerSpectralDensity, time::AbstractVector{<:Real}, name::Symbol
)
    Δω = psd.ω[2] - psd.ω[1]

    A = sqrt.(2 * psd.p * Δω)
    A[iszero.(psd.ω)] .= 0.0

    return SpectralRepresentation(
        psd,
        time,
        psd.ω * time',
        Δω,
        A,
        name,
        [Symbol("$(name)_$i") for i in 1:length(psd.ω)],
    )
end

function sample(sr::SpectralRepresentation, n::Integer=1)
    return DataFrame(names(sr) .=> eachcol(rand(Uniform(0, 2π), (n, length(sr.psd.ω)))))
end

function evaluate(sr::SpectralRepresentation, ϕ::AbstractVector{<:Real})
    return sqrt(2) * vec(sr.A' * cos.(sr.ωt .+ ϕ))
end

function (sr::SpectralRepresentation)(ϕ::AbstractVector{<:Real})
    return evaluate(sr, ϕ)
end

# function evaluate(sr::SpectralRepresentation, row::DataFrameRow)
#     return evaluate(sr, collect(row[sr.ϕnames]))
# end

# function evaluate(sr::SpectralRepresentation, row::DataFrameRow)
#     return evaluate(sr, collect(row[sr.ϕnames]))
# end

# function (sr::SpectralRepresentation)(row::DataFrameRow)
#     return evaluate(sr, row)
# end

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
    return sr.ϕnames
end
