"""
    GaussianMixtureModel(names::Vector{Symbol}, number_components::Integer)

    Creates a Gaussian Mixture Model (GMM) with the specified number of components and dimensions.
    The GMM is initialized with multivariate normal distributions for each component, with means set to zero and identity covariance matrices.
"""
struct GaussianMixtureModel <:RandomUQInput
    mixture::MixtureModel
    names::Vector{Symbol}

    function GaussianMixtureModel(
        data::DataFrame,
        number_components::Integer;
        maximum_iterations::Integer=100,
        tolerance::Number=1e-4
    )
        # Fit mixture model to the data using EM algorithm
        mixture = fit_gaussian_mixture(number_components, Matrix(data);
            maximum_iterations=maximum_iterations, tolerance=tolerance)

        return new(mixture, Symbol.(names(data)))
    end

    function GaussianMixtureModel(mixture::MixtureModel, names::Vector{Symbol})
        if length(mixture.components[1].μ) != length(names)
            throw(ArgumentError("Mismatch between dimensions and names."))
        end
        return new(mixture, names)
    end
end

"""
    sample(gmm::GaussianMixtureModel, n::Integer=1)

Generates `n` samples from the Gaussian Mixture Model. Returns a DataFrame.
"""
function sample(gmm::GaussianMixtureModel, n::Integer=1)
    samples = rand(gmm.mixture, n)
    # Handle one-dimensional case
    if length(gmm.names) == 1
        # If the GMM has only one dimension, return a DataFrame with a single column
        return DataFrame(gmm.names[1] => samples)
    end
    return DataFrame(transpose(samples), gmm.names)
end

function to_standard_normal_space!(gmm::GaussianMixtureModel, df::DataFrame)
	error("GaussianMixtureModels cannot be mapped to standard normal space.")
end

names(gmm::GaussianMixtureModel) = gmm.names
dimensions(gmm::GaussianMixtureModel) = length(gmm.names)

mean(gmm::GaussianMixtureModel) = mean(gmm.mixture)
var(gmm::GaussianMixtureModel) = var(gmm.mixture)
pdf(gmm::GaussianMixtureModel, x::Union{Vector{<:Real}, <:Real}) = pdf(gmm.mixture, x)
logpdf(gmm::GaussianMixtureModel, x::Union{Vector{<:Real}, <:Real}) = logpdf(gmm.mixture, x)
minimum(gmm::GaussianMixtureModel) = minimum(gmm.mixture)
maximum(gmm::GaussianMixtureModel) = maximum(gmm.mixture)
insupport(gmm::GaussianMixtureModel, x::Union{Vector{<:Real}, <:Real}) = insupport(gmm.mixture, x)

"""
    fit_gaussian_mixture(number_components::Integer, number_samples::Integer, data::Matrix; maximum_iterations::Integer=100, tolerance::Number=1e-4)

Fits a Gaussian Mixture Model (GMM) to the given data using the Expectation-Maximization (EM) algorithm.

# Arguments
- `number_components::Integer`: Number of mixture components (clusters).
- `number_samples::Integer`: Number of data points (rows in `X`).
- `data::Matrix`: Data matrix of size `N × D`, where each row is a data point and each column is a feature.
- `maximum_iterations::Integer=100`: Maximum number of EM iterations. Defaults to 100.
- `tolerance::Number=1e-4`: (Optional) Convergence tolerance for the change in log-likelihood. Defaults to 1e-4.

# Returns
- A `Distributions.MixtureModel` object containing the fitted Gaussian components and their mixing weights.
"""
function fit_gaussian_mixture(number_components::Integer, data::Matrix; maximum_iterations::Integer=100, tolerance::Number=1e-4)
    # Check inputs
    if number_components <= 0
        throw(ArgumentError("Number of components must be a positive integer."))
    end
    if size(data, 2) < 2
        throw(ArgumentError("Input must be at least two-dimensional to fit GMM."))
    end

    # Initialize parameters
    number_samples, _ = size(data)
    μ, Σ, π = _initialize_gaussian_mixture_model(data, number_components)
    log_likelihood_old = -Inf

    for _ in 1:maximum_iterations
        γ = _expectation_step(data, μ, Σ, π)
        μ, Σ, π = _maximization_step(data, γ)

        # Log-likelihood (for convergence check)
        log_likelihood = sum(log(sum(π[k] * pdf(MvNormal(μ[k, :], Σ[k]), data[i, :]) for k in 1:number_components)) for i in 1:number_samples)

        abs(log_likelihood - log_likelihood_old) < tolerance ? break : continue
        log_likelihood_old = log_likelihood
    end

    return MixtureModel([MvNormal(μ[k, :], Σ[k]) for k in 1:number_components], π)

end

function _initialize_gaussian_mixture_model(data::Matrix{Float64}, number_components::Integer)
    # Randomly initialize parameters for GMM
    number_samples, dimensions = size(data)
    μ = data[rand(1:number_samples, number_components), :]  # random initialization of means
    Σ = [I(dimensions) for _ in 1:number_components]        # identity covariance matrices
    π = fill(1/number_components, number_components)        # equal mixing weights
    return μ, Σ, π
end

function _expectation_step(data::Matrix{Float64}, μ::Matrix{Float64}, Σ::Vector, π::Vector{Float64})
    number_samples, number_components = size(data, 1), length(π)
    # Compute responsibilities and normalize
    γ = [π[k] * pdf(MvNormal(μ[k, :], Σ[k]), data[i, :]) for i in 1:number_samples, k in 1:number_components]
    γ = γ ./ sum(γ, dims=2)
    return γ
end

function _maximization_step(data::Matrix{Float64}, γ::Matrix{Float64})
    number_samples, dimensions = size(data)
    number_components = size(γ, 2)
    Nₖ = sum(γ, dims=1)             # effective number of points in cluster 𝑘
    π = vec(Nₖ ./ number_samples)   # weights: relative number of points in each cluster
    μ = γ' * data ./ Nₖ'            # means: weighted average of points in each cluster
    Σ = Vector{Matrix{Float64}}(undef, number_components)
    for k in 1:number_components
        data_centered = data .- μ[k:k, :]
        Σₖ = zeros(dimensions, dimensions)
        for i in 1:number_samples
            Σₖ += γ[i, k] * (data_centered[i, :] * data_centered[i, :]')
        end
        Σ[k] = Σₖ / Nₖ[k]
    end
    return μ, Σ, π
end
