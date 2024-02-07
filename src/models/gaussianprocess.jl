"""
    GaussianProcess(data::DataFrame, dependendVarName::Symbol, deg::Int, dim::Int)

Creates a gaussian process prior ....

# Examples
```jldoctest
julia> data = DataFrame(x = 1:10, y = [1, 4, 10, 15, 24, 37, 50, 62, 80, 101]);

julia> rs = ResponseSurface(data, :y, 2) |> DisplayAs.withcontext(:compact => true)
ResponseSurface([0.483333, -0.238636, 1.01894], :y, [:x], 2, Monomial{Commutative{CreationOrder}, Graded{LexOrder}}[1, x₁, x₁²])
```
"""
mutable struct GaussianProcessRegressor <: UQModel
    gp::GPBase
    inputs::Union{Vector{<:UQInput}, Vector{Symbol}}
    output::Symbol
    input_normalizer::Union{Normalizer, Nothing}
    output_normalizer::Union{Normalizer, Nothing}
end

struct Normalizer
    μ::Union{Real, Matrix{<:Real}}
    σ::Union{Real, Matrix{<:Real}}
end

normalize(data::Union{Vector{Real}, Matrix{<:Real}}, normalizer::Normalizer) = (data .- normalizer.μ) ./ normalizer.σ
denormalize(data::Union{Vector{Real}, Matrix{<:Real}}, normalizer::Normalizer) = data .* normalizer.σ .+ normalizer.μ

function gaussianprocess(
    df::DataFrame,
    inputs::Vector{Symbol},
    output::Symbol,
    kernel::Kernel,
    mean::Mean=MeanZero(),
    log_noise::Real=-2.0,
    normalize_input::Bool=false,
    normalize_output::Bool=false
)
    X = Matrix(df[:, inputs])'
    y = df[:, output]

    if normalize_input
        normalizer_in = Normalizer(mean(X, dims=2), std(X, dims=2))
        X = normalize(X, normalizer_in)
    else
        normalizer_in = nothing
    end
    if normalize_output
        normalizer_out = Normalizer(mean(y), std(y))
        y = normalize(y, normalizer_out)
    else
        normalizer_out = nothing
    end
    
    gp = GP(X, y, mean, kernel)
    optimize!(gp)

    gp = GaussianProcess(
        gp, inputs, output, size(X, 2), InputStandardizationGP(minimum(X), maximum(X))
    )

    return gp, df
end



