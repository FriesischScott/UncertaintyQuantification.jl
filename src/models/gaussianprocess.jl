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
    method::Union{Optim.LBFGS, Optim.ConjugateGradient} # not sure how or even if to support multiple solvers
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

function polynomialchaos(
    inputs::Vector{<:UQInput},
    model::Vector{<:UQModel},
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    _::GaussQuadrature,
)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)
    random_names = names(random_inputs)

    nodes = mapreduce(
        n -> [n...]', vcat, Iterators.product(quadrature_nodes.(Ψ.p + 1, Ψ.bases)...)
    )
    weights = map(prod, Iterators.product(quadrature_weights.(Ψ.p + 1, Ψ.bases)...))

    samples = DataFrame(map_from_bases(Ψ, nodes), random_names)
    to_physical_space!(random_inputs, samples)

    if !isempty(deterministic_inputs)
        samples = hcat(samples, sample(deterministic_inputs, size(nodes, 1)))
    end

    evaluate!(model, samples)

    y = mapreduce(
        (x, w, f) -> f * w * evaluate(Ψ, collect(x)),
        +,
        eachrow(nodes),
        weights,
        samples[:, output],
    )

    return PolynomialChaosExpansion(y, Ψ, output, random_inputs), samples
end

function polynomialchaos(
    inputs::UQInput,
    model::Vector{<:UQModel},
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    gq::GaussQuadrature,
)
    return polynomialchaos([inputs], model, Ψ, output, gq)
end

function polynomialchaos(
    inputs::Vector{<:UQInput},
    model::UQModel,
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    gq::GaussQuadrature,
)
    return polynomialchaos(inputs, [model], Ψ, output, gq)
end

function polynomialchaos(
    inputs::UQInput,
    model::UQModel,
    Ψ::PolynomialChaosBasis,
    output::Symbol,
    gq::GaussQuadrature,
)
    return polynomialchaos([inputs], [model], Ψ, output, gq)
end

# what should this return?
function evaluate!(gpr::GaussianProcessRegressor, df::DataFrame) # this now gives mean and variance at inputs
    data = Matrix(df[:, names(gpr.inputs)])'
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

    samps = DataFrame(map_from_bases(pce.Ψ, samps), names(pce.inputs))
    to_physical_space!(pce.inputs, samps)

    samps[!, pce.output] = out
    return samps
end

