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
struct Normalizer # Maybe move this somewhere else...
    μ::Union{Real, Matrix{<:Real}}
    σ::Union{Real, Matrix{<:Real}}
end

normalize(data::Union{Vector{Real}, Matrix{<:Real}}, normalizer::Normalizer) = (data .- normalizer.μ) ./ normalizer.σ
denormalize(data::Union{Vector{Real}, Matrix{<:Real}}, normalizer::Normalizer) = data .* normalizer.σ .+ normalizer.μ

mutable struct GaussianProcessRegressor <: UQModel
    gp::GPBase
    inputs::Union{Vector{<:UQInput}, Vector{Symbol}}
    output::Symbol
    input_normalizer::Union{Normalizer, Nothing}
    output_normalizer::Union{Normalizer, Nothing}
end

function gaussianprocess(
    df::DataFrame,
    inputs::Vector{Symbol},
    output::Symbol,
    kernel::Kernel,
    mean::GaussianProcesses.Mean=MeanZero(),
    log_noise::Real=-2.0,
    normalize_input::Bool=false,
    normalize_output::Bool=false
)
    X = Matrix(df[:, inputs])'
    y = df[:, output]

    if normalize_input
        input_normalizer = Normalizer(mean(X, dims=2), std(X, dims=2))
        X = normalize(X, input_normalizer)
    else
        input_normalizer = nothing
    end
    if normalize_output
        output_normalizer = Normalizer(mean(y), std(y))
        y = normalize(y, output_normalizer)
        log_noise -= log(output_normalizer.σ)
    else
        output_normalizer = nothing
    end
    
    gp = GP(X, y, mean, kernel, log_noise)

    gp = GaussianProcessRegressor(
        gp, inputs, output, 
        input_normalizer, output_normalizer
        )

    return gp, df
end

# Wrapper for optimize! method from GaussianProcesses.jl
function optimize_hyperparams!(gpr::GaussianProcessRegressor, args...; method = LBFGS(), 
    domean::Bool = true, kern::Bool = true, noise::Bool = true, 
    lik::Bool = true, meanbounds = nothing, kernbounds = nothing,
    noisebounds = nothing, likbounds = nothing, kwargs...
)

    optimize!(gpr.gp, args...; method = method, 
    domean=domean, kern=kern, noise=noise, lik=lik,
    meanbounds=meanbounds, kernbounds=kernbounds,
    noisebounds=noisebounds, likbounds=likbounds, 
    kwargs...)
end



