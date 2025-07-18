struct GaussianProcess <: UQModel
    gp::AbstractGPs.PosteriorGP
    input::Union{Symbol, Vector{Symbol}}
    output::Symbol
    datatransformer::DataTransform
end

# Build from Dataframe
function GaussianProcess(
    gp::Union{AbstractGPs.GP, NoisyGP},
    data::DataFrame,
    output::Symbol,
    instandard::StandardizeInput,
    outstandard::StandardizeOutput,
    optimization::AbstractHyperparameterOptimization
) # should we use keyword args?
    input = propertynames(data[:, Not(output)])
    datatransformer = DataTransform(
        data, instandard, input, outstandard, output
    )
    x = datatransformer.input(data)
    y = datatransformer.output(data)
    optimized_gp = optimize_hyperparameters(gp, x, y, optimization)
    posterior_gp = posterior(optimized_gp(x), y)
    return GaussianProcess(
        posterior_gp,
        input,
        output,
        datatransformer
    )
end

function GaussianProcess(
    gp::AbstractGPs.GP,
    data::DataFrame,
    output::Symbol,
    instandard::StandardizeInput,
    optimization::AbstractHyperparameterOptimization
) # should we use keyword args?
    return GaussianProcess(
        gp,
        data,
        output,
        instandard,
        StandardizeOutput(false),
        optimization
    )
end

function GaussianProcess(
    gp::AbstractGPs.GP,
    data::DataFrame,
    output::Symbol,
    outstandard::StandardizeOutput,
    optimization::AbstractHyperparameterOptimization
) # should we use keyword args?
    return GaussianProcess(
        gp,
        data,
        output,
        StandardizeInput(false),
        outstandard,
        optimization
    )
end

function GaussianProcess(
    gp::AbstractGPs.GP,
    data::DataFrame,
    output::Symbol,
    optimization::AbstractHyperparameterOptimization
) # should we use keyword args?
    return GaussianProcess(
        gp,
        data,
        output,
        StandardizeInput(false),
        StandardizeOutput(false),
        optimization
    )
end

function GaussianProcess(
    gp::AbstractGPs.GP,
    data::DataFrame,
    output::Symbol
) # should we use keyword args?
    return GaussianProcess(
        gp,
        data,
        output,
        StandardizeInput(false),
        StandardizeOutput(false),
        NoOptimization()
    )
end

# Build with UQmodel
function GaussianProcess(
    gp::Union{AbstractGPs.GP, NoisyGP},
    input::Union{UQInput, Vector{<:UQInput}},
    model::Union{UQModel, Vector{<:UQModel}},
    output::Symbol,
    experimentaldesign::Union{AbstractMonteCarlo, AbstractDesignOfExperiments},
    instandard::StandardizeInput,
    outstandard::StandardizeOutput,
    optimization::AbstractHyperparameterOptimization
)
    data = sample(input, experimentaldesign) # need to be able to pass experimental design
    evaluate!(model, data)

    datatransformer = DataTransform(
        data, instandard, input, outstandard, output
    )
    x = datatransformer.input(data)
    y = datatransformer.output(data)
    optimized_gp = optimize_hyperparameters(gp, x, y, optimization)
    posterior_gp = posterior(optimized_gp(x), y)
    return GaussianProcess(
        posterior_gp,
        names(input),
        output,
        datatransformer
    )
end

function GaussianProcess(
    gp::Union{AbstractGPs.GP, NoisyGP},
    input::Union{UQInput, Vector{<:UQInput}},
    model::Union{UQModel, Vector{<:UQModel}},
    output::Symbol,
    experimentaldesign::Union{AbstractMonteCarlo, AbstractDesignOfExperiments},
    instandard::StandardizeInput,
    optimization::AbstractHyperparameterOptimization
)
    return GaussianProcess(
        gp,
        input,
        model, 
        output,
        experimentaldesign,
        instandard,
        StandardizeOutput(false),
        optimization
    )
end

function GaussianProcess(
    gp::Union{AbstractGPs.GP, NoisyGP},
    input::Union{UQInput, Vector{<:UQInput}},
    model::Union{UQModel, Vector{<:UQModel}},
    output::Symbol,
    experimentaldesign::Union{AbstractMonteCarlo, AbstractDesignOfExperiments},
    outstandard::StandardizeOutput,
    optimization::AbstractHyperparameterOptimization
)
    return GaussianProcess(
        gp,
        input,
        model, 
        output,
        experimentaldesign,
        StandardizeInput(false),
        outstandard,
        optimization
    )
end

function GaussianProcess(
    gp::Union{AbstractGPs.GP, NoisyGP},
    input::Union{UQInput, Vector{<:UQInput}},
    model::Union{UQModel, Vector{<:UQModel}},
    output::Symbol,
    experimentaldesign::Union{AbstractMonteCarlo, AbstractDesignOfExperiments},
    optimization::AbstractHyperparameterOptimization
)
    return GaussianProcess(
        gp,
        input,
        model, 
        output,
        experimentaldesign,
        StandardizeInput(false),
        StandardizeOutput(false),
        optimization
    )
end

function GaussianProcess(
    gp::Union{AbstractGPs.GP, NoisyGP},
    input::Union{UQInput, Vector{<:UQInput}},
    model::Union{UQModel, Vector{<:UQModel}},
    output::Symbol,
    experimentaldesign::Union{AbstractMonteCarlo, AbstractDesignOfExperiments}
)
    return GaussianProcess(
        gp,
        input,
        model, 
        output,
        experimentaldesign,
        StandardizeInput(false),
        StandardizeOutput(false),
        NoOptimization()
    )
end

# what should this calculate? Calculates only mean for now
function evaluate!(gp::GaussianProcess, data::DataFrame)
    x = gp.datatransformer.input(data)
    y = mean(gp.gp(x))

    data[!, gp.output] = inverse_transform(y, gp.datatransformer.output)
    return nothing
end