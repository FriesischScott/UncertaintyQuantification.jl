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
    input::Union{Vector{<:UQInput}, Vector{Symbol}}
    output::Symbol
    input_normalizer::Union{ZScoreTransform, Nothing}
    output_normalizer::Union{ZScoreTransform, Nothing}
end

function normalize!(
    input::Union{Vector{<:Real}, Matrix{<:Real}},
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

struct Optimizer
    # maybe give one default and allow JuMP structs
    method::Union{Optim.LBFGS, Optim.ConjugateGradient} # not sure how or even if to support multiple solvers (Matlab uses QuasiNewton default)
    optim_options::Dict # maybe there is a better option than using dicts for this
    hyperparams::Dict
    bounds::Dict
end

struct ExperimentalDesign # not sure about the name
    sim::AbstractMonteCarlo # could also allow doe
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
    optimizer::Union{Optimizer, Nothing}=Optimizer(),
    normalize_input::Bool=false,
    normalize_output::Bool=false
)
    x = copy(Matrix(df[:, input])')
    y = df[:, output]
    input_normalizer, output_normalizer, log_noise = normalize!(x, y, normalize_input, normalize_output, log_noise)
    
    gp = GP(x, y, mean, kernel, log_noise)
    if !isnothing(optimizer)
        optimize!(gp; 
            method=optimizer.method, 
            optimizer.hyperparams...,
            optimizer.bounds...,
            optimizer.optim_options...
        )
    end

    gp = GaussianProcessRegressor(
        gp, input, output, 
        input_normalizer, output_normalizer
        )

    return gp, df
end

function gaussianprocess(
    input::Vector{<:UQInput},
    model::Vector{<:UQModel},
    output::Symbol,
    ed::ExperimentalDesign,
    kernel::Kernel,
    mean::GaussianProcesses.Mean=MeanZero(),
    log_noise::Real=-2.0,
    optimizer::Union{Optimizer, Nothing}=Optimizer(),
    normalize_output::Bool=false
)
    samples = sample(input, ed.sim)
    evaluate!(model, samples)

    random_input = filter(i -> isa(i, RandomUQInput), input)
    random_names = names(random_input)

    to_standard_normal_space!(random_input, samples) # not sure if this is save to do in every case
    x = copy(Matrix(samples[:, random_names])')
    y = df[:, output]
    _, output_normalizer, log_noise = normalize!(x, y, false, normalize_output, log_noise) # do not need input normalizer here
    
    gp = GP(x, y, mean, kernel, log_noise)
    if !isnothing(optimizer)
        optimize!(gp; 
            method=optimizer.method, 
            optimizer.hyperparams...,
            optimizer.bounds...,
            optimizer.optim_options...
        )
    end

    gp = GaussianProcessRegressor(
        gp, input, output, 
        _, output_normalizer
        )
    to_physical_space!(random_input, samples)

    return gp, samples
end

function gaussianprocess(
    input::UQInput,
    model::Vector{<:UQModel},
    output::Symbol,
    ed::ExperimentalDesign,
    kernel::Kernel,
    mean::GaussianProcesses.Mean=MeanZero(),
    log_noise::Real=-2.0,
    optimizer::Union{Optimizer, Nothing}=Optimizer(),
    normalize_output::Bool=false
)
    return gaussianprocess(
        [input], model, output, 
        ed, kernel, mean, log_noise, 
        optimizer, normalize_output
    )
end

function gaussianprocess(
    input::Vector{<:UQInput},
    model::UQModel,
    output::Symbol,
    ed::ExperimentalDesign,
    kernel::Kernel,
    mean::GaussianProcesses.Mean=MeanZero(),
    log_noise::Real=-2.0,
    optimizer::Union{Optimizer, Nothing}=Optimizer(),
    normalize_output::Bool=false
)
    return gaussianprocess(
        input, [model], output, 
        ed, kernel, mean, log_noise, 
        optimizer, normalize_output
    )
end

function gaussianprocess(
    input::UQInput,
    model::UQModel,
    output::Symbol,
    ed::ExperimentalDesign,
    kernel::Kernel,
    mean::GaussianProcesses.Mean=MeanZero(),
    log_noise::Real=-2.0,
    optimizer::Union{Optimizer, Nothing}=Optimizer(),
    normalize_output::Bool=false
)
    return gaussianprocess(
        [input], [model], output, 
        ed, kernel, mean, log_noise, 
        optimizer, normalize_output
    )
end

# what should this return?
function evaluate!(gpr::GaussianProcessRegressor, df::DataFrame) # this now gives mean and variance at input
    data = Matrix(df[:, names(gpr.input)])'
    if !isnothing(gpr.input_normalizer)
        μ, Σ = predict_y(gpr.gp, StatsBase.transform!(grp.input_normalizer, data))
    else
        μ, Σ = predict_y(gpr.gp, data)
    end

    if !isnothing(grp.output_normalizer)
        μ[:] = μ .* gpr.output_normalizer.scale[1] .+ gpr.output_normalizer.mean[1] 
        Σ[:] = Σ .* gpr.output_normalizer.scale[1]^2
    end

    df[!, Symbol(gpr.output, "_mean")] = μ
    df[!, Symbol(gpr.output, "_var")] = Σ
    return nothing
end

function sample(pce::PolynomialChaosExpansion, n::Integer)
    samps = hcat(sample.(n, pce.Ψ.bases)...)
    out = map(row -> dot(pce.y, evaluate(pce.Ψ, collect(row))), eachrow(samps))

    samps = DataFrame(map_from_bases(pce.Ψ, samps), names(pce.input))
    to_physical_space!(pce.input, samps)

    samps[!, pce.output] = out
    return samps
end

