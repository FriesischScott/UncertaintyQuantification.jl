struct GaussianProcess <: UQModel
    gp::AbstractGPs.PosteriorGP
    input::Union{Symbol, Vector{Symbol}}
    output::Symbol
    standardizer::DataStandardizer
end

# ---------------- Build from DataFrame ----------------
function GaussianProcess(
    gp::Union{AbstractGPs.GP, NoisyGP},
    data::DataFrame,
    output::Symbol;
    input_transform::AbstractDataTransform=IdentityTransform(),
    output_transform::AbstractDataTransform=IdentityTransform(),
    optimization::AbstractHyperparameterOptimization=MaximumLikelihoodEstimation()
) 
    input = propertynames(data[:, Not(output)]) # Is this always the case?

    # build in- and output transforms
    dts = DataStandardizer(
        data, input, output, 
        InputTransform(input_transform), 
        OutputTransform(output_transform)
    )

    # transform data
    x = dts.fᵢ(data)
    y = dts.fₒ(data)

    # build posterior gp
    optimized_gp = optimize_hyperparameters(gp, x, y, optimization)
    posterior_gp = posterior(optimized_gp(x), y)
    return GaussianProcess(
        posterior_gp,
        input,
        output,
        dts
    )
end

# ---------------- Build from UQModel ----------------
function GaussianProcess(
    gp::Union{AbstractGPs.GP, NoisyGP},
    input::Union{UQInput, Vector{<:UQInput}},
    model::Union{UQModel, Vector{<:UQModel}},
    output::Symbol,
    experimentaldesign::Union{AbstractMonteCarlo, AbstractDesignOfExperiments};
    input_transform::AbstractDataTransform=IdentityTransform(),
    output_transform::AbstractDataTransform=IdentityTransform(),
    optimization::AbstractHyperparameterOptimization=MaximumLikelihoodEstimation()
)
    # build DataFrame
    data = sample(input, experimentaldesign)
    evaluate!(model, data)

    # build in- and output transforms
    dts = DataStandardizer(
        data, input, output, 
        InputTransform(input_transform), 
        OutputTransform(output_transform)
    )

    # transform data
    x = dts.fᵢ(data)
    y = dts.fₒ(data)

    # build posterior gp
    optimized_gp = optimize_hyperparameters(gp, x, y, optimization)
    posterior_gp = posterior(optimized_gp(x), y)

    return GaussianProcess(
        posterior_gp,
        names(input),
        output,
        dts
    )
end

# what should this calculate? Calculates only mean for now
function evaluate!(gp::GaussianProcess, data::DataFrame)
    x = gp.standardizer.fᵢ(data)
    y = mean(gp.gp(x))

    data[!, gp.output] = gp.standardizer.fₒ⁻¹(y) # applying inverse transform to output
    return nothing
end