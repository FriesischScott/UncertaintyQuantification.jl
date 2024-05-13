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
    input_normalizer::Union{ZScoreTransform, Nothing}
    output_normalizer::Union{ZScoreTransform, Nothing}
end

function normalize!( # maybe different name as this is only supposed to be used in the GP context
    input::Union{Vector{<:Real}, AbstractMatrix{<:Real}},
    output::Vector{<:Real},
    normalize_input::Bool,
    normalize_output::Bool,
    log_noise::Real
)
    if normalize_input
        input_normalizer = fit(ZScoreTransform, input)
        input[:] = StatsBase.transform(input_normalizer, input)
    else
        input_normalizer = nothing
    end

    if normalize_output
        output_normalizer = fit(ZScoreTransform, output)
        output[:] = StatsBase.transform(output_normalizer, output)
        log_noise -= log(output_normalizer.scale[1])
    else
        output_normalizer = nothing
    end

    return input_normalizer, output_normalizer, log_noise
end

struct Optimizer # there is probably a better way to design this
    optimizer
    opt_kwargs::Dict
    hyperparams::Dict
    bounds::Dict
    # should I add number of optimizer runs?
end

Optimizer() = Optimizer(
            LBFGS(), 
            Dict(), 
            Dict(:domean => true, :kern => true, :noise => true, :lik => true),
            Dict(:meanbounds => nothing, :kernbounds => nothing, :noisebounds => nothing, :likbounds => nothing)
            )

function gaussianprocess(
    df::DataFrame,
    input::Vector{Symbol},
    output::Symbol,
    kernel::Kernel,
    mean::GaussianProcesses.Mean=MeanZero(),
    log_noise::Real=-2.0,
    optimizer::Union{Optimizer, Nothing}=Optimizer(), # there is probably a better way to design this
    normalize_input::Bool=false,
    normalize_output::Bool=false
)
    x = Matrix(df[:, input])'
    y = df[:, output]

    input_normalizer, output_normalizer, log_noise = normalize!(x, y, normalize_input, normalize_output, log_noise)
    
    gp = GP(x, y, mean, kernel, log_noise)

    if !isnothing(optimizer)
        optimize!(gp; 
        method=optimizer.optimizer,
        optimizer.hyperparams...,
        optimizer.bounds...,
        optimizer.opt_kwargs... 
        )
    end

    gp = GaussianProcessRegressor(
        gp, input, output, 
        input_normalizer, output_normalizer
        )

    return gp, df # this method does not really need to return df
end

# Wrapper for optimize! method from GaussianProcesses.jl
# function optimize_hyperparams!(gpr::GaussianProcessRegressor, args...; method = LBFGS(), 
#     domean::Bool = true, kern::Bool = true, noise::Bool = true, 
#     lik::Bool = true, meanbounds = nothing, kernbounds = nothing,
#     noisebounds = nothing, likbounds = nothing, kwargs...
# )

#     optimize!(gpr.gp, args...; method = method, 
#     domean=domean, kern=kern, noise=noise, lik=lik,
#     meanbounds=meanbounds, kernbounds=kernbounds,
#     noisebounds=noisebounds, likbounds=likbounds, 
#     kwargs...)
# end



