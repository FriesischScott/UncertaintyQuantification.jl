"""
    GaussianMixtureModel(names::Vector{Symbol}, number_components::Integer)

    Creates a Gaussian Mixture Model (GMM) with the specified number of components and dimensions.
    The GMM is initialized with multivariate normal distributions for each component, with means set to zero and identity covariance matrices.
"""
mutable struct GaussianMixtureModel <:RandomUQInput
    mixture::MixtureModel
    names::Vector{Symbol}

    function GaussianMixtureModel(names::Vector{Symbol}, number_components::Integer)
        dimension = length(names)
        mixture = MixtureModel(
            [MvNormal(zeros(dimension), I(dimension)) for k in 1:number_components], fill(1/number_components, number_components)
        )
        return new(mixture, names)
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
    return DataFrame(transpose(samples), gmm.names)
end

function to_standard_normal_space!(gmm::GaussianMixtureModel, df::DataFrame)
	error("GaussianMixtureModels cannot be mapped to standard normal space.")
end

names(gmm::GaussianMixtureModel) = gmm.names
dimensions(gmm::GaussianMixtureModel) = length(gmm.names)

mean(gmm::GaussianMixtureModel) = mean(gmm.mixture)
var(gmm::GaussianMixtureModel) = var(gmm.mixture)
pdf(gmm::GaussianMixtureModel, x::Vector{Float64}) = pdf(gmm.mixture, x)
logpdf(gmm::GaussianMixtureModel, x::Vector{Float64}) = logpdf(gmm.mixture, x)
minimum(gmm::GaussianMixtureModel) = minimum(gmm.mixture)
maximum(gmm::GaussianMixtureModel) = maximum(gmm.mixture)
insupport(gmm::GaussianMixtureModel, x::Vector{Float64}) = insupport(gmm.mixture, x)

"""
    fit!(gmm::GaussianMixtureModel, data::DataFrame; max_iter::Integer=100, tol::Number=1e-4)

Fits a Gaussian Mixture Model to the given data using the Expectation-Maximization algorithm.

# Arguments
- `gmm::GaussianMixtureModel`: The GMM object to be fitted. The model parameters will be updated in-place.
- `data::DataFrame`: The input data as a DataFrame, where each row corresponds to a data point.
- `max_iter::Integer=100`: (Optional) Maximum number of EM iterations. Defaults to 100.
- `tol::Number=1e-4`: (Optional) Convergence tolerance for the change in log-likelihood. Defaults to 1e-4.

# Example
```julia
df = DataFrame(x1=randn(100), x2=randn(100))
gmm = GaussianMixtureModel([:x1, :x2], 3)
fit!(gmm, df)
```
"""
function fit!(gmm::GaussianMixtureModel, data::DataFrame; max_iter::Integer=100, tol::Number=1e-4)
    # Initialize parameters
    if ncol(data) != length(gmm.names)
        throw(ArgumentError("Number of columns in data must match the number of names in the GMM."))
    end
    X = Matrix(data)
    K = length(gmm.mixture.components)
    N = nrow(data)
    μ, Σ, π = _initialize_gaussian_mixture_model(X, K)

    log_likelihood_old = -Inf

    for _ in 1:max_iter
        γ = _expectation_step(X, μ, Σ, π)
        μ, Σ, π = _maximization_step(X, γ)

        # Log-likelihood (for convergence check)
        log_likelihood = sum(log(sum(π[k] * pdf(MvNormal(μ[k, :], Σ[k]), X[i, :]) for k in 1:K)) for i in 1:N)

        abs(log_likelihood - log_likelihood_old) < tol ? break : continue
        log_likelihood_old = log_likelihood
    end

    gmm.mixture = MixtureModel([MvNormal(μ[k, :], Σ[k]) for k in 1:K], π)

end

function _initialize_gaussian_mixture_model(X::Matrix{Float64}, K::Integer)
    # Randomly initialize parameters for GMM
    N, D = size(X)
    μ = X[rand(1:N, K), :]      # random initialization of means
    Σ = [I(D) for _ in 1:K]     # identity covariance matrices
    π = fill(1/K, K)            # equal mixing weights
    return μ, Σ, π
end

function _expectation_step(X::Matrix{Float64}, μ::Matrix{Float64}, Σ::Vector, π::Vector{Float64})
    N, K = size(X, 1), length(π)
    # Compute responsibilities and normalize
    γ = [π[k] * pdf(MvNormal(μ[k, :], Σ[k]), X[i, :]) for i in 1:N, k in 1:K]
    γ = γ ./ sum(γ, dims=2)
    return γ
end

function _maximization_step(X::Matrix{Float64}, γ::Matrix{Float64})
    N, D = size(X)
    K = size(γ, 2)
    Nₖ = sum(γ, dims=1) # effective number of points in cluster 𝑘
    π = vec(Nₖ ./ N)    # weights: relative number of points in each cluster
    μ = γ' * X ./ Nₖ'   # means: weighted average of points in each cluster
    Σ = Vector{Matrix{Float64}}(undef, K)
    for k in 1:K
        X_centered = X .- μ[k:k, :]
        Σₖ = zeros(D, D)
        for i in 1:N
            Σₖ += γ[i, k] * (X_centered[i, :] * X_centered[i, :]')
        end
        Σ[k] = Σₖ / Nₖ[k]
    end
    return μ, Σ, π
end
