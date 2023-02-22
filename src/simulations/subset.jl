abstract type AbstractSubSetSimulation end

"""
    SubSetSimulation(n::Integer, target::Float64, levels::Integer, proposal::UnivariateDistribution)

Defines the properties of a Subset simulation where `n` is the number of initial samples,
`target` is the target probability of failure at each level, `levels` is the maximum number
of levels and `proposal` is the proposal distribution for the markov chain monte carlo.

# Examples

```jldoctest
julia> SubSetSimulation(100, 0.1, 10, Uniform(-0.2, 0.2))
SubSetSimulation(100, 0.1, 10, Uniform{Float64}(a=-0.2, b=0.2))
```

# References

[auEstimationSmallFailure2001](@cite)
"""
struct SubSetSimulation <: AbstractSubSetSimulation
    n::Integer
    target::Float64
    levels::Integer
    proposal::UnivariateDistribution

    function SubSetSimulation(
        n::Integer, target::Float64, levels::Integer, proposal::UnivariateDistribution
    )
        skewness(proposal) != 0.0 && error("proposal must be a symmetric distribution")
        mean(proposal) != median(proposal) &&
            error("proposal must be a symmetric distribution")
        return new(n, target, levels, proposal)
    end
end

"""
    SubSetInfinity(n::Integer, target::Float64, levels::Integer, s::Real)

Defines the properties of a Subset-∞ simulation where `n` is the number of initial samples,
`target` is the target probability of failure at each level, `levels` is the maximum number
of levels and `s` is the standard deviation for the proposal samples.

# Examples

```jldoctest
julia> SubSetInfinity(100, 0.1, 10, 0.5)
SubSetInfinity(100, 0.1, 10, 0.5)
```

# References

[auRareEventSimulation2016](@cite)

[patelliEfficientMonteCarlo2015](@cite)
"""
struct SubSetInfinity <: AbstractSubSetSimulation
    n::Integer
    target::Float64
    levels::Integer
    s::Real

    function SubSetInfinity(n::Integer, target::Float64, levels::Integer, s::Real)
        (0 <= s <= 1) || error("standard deviation must be between 0.0")
        return new(n, target, levels, s)
    end
end

function sample(inputs::Array{<:UQInput}, sim::AbstractSubSetSimulation)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)

    n_rv = count_rvs(random_inputs)
    rv_names = names(random_inputs)

    samples = DataFrame(rv_names .=> eachcol(rand(Normal(), sim.n, n_rv)))

    to_physical_space!(random_inputs, samples)

    if !isempty(deterministic_inputs)
        samples = hcat(samples, sample(deterministic_inputs, sim.n))
    end
    return samples
end

function probability_of_failure(
    models::Union{Array{<:UQModel},UQModel},
    performancefunction::Function,
    inputs::Union{Array{T},T} where {T<:UQInput},
    sim::AbstractSubSetSimulation,
)
    samples = [sample(inputs, sim)]

    evaluate!(models, samples[end])

    performance = [performancefunction(samples[end])]

    number_of_seeds = Int64(max(1, ceil(sim.n * sim.target)))
    samples_per_seed = Int64(floor(sim.n / number_of_seeds))

    threshold = zeros(sim.levels, 1)
    pf = ones(sim.levels, 1)
    cov = zeros(sim.levels, 1)

    for i in 1:(sim.levels)
        sorted_performance = sort(performance[end])
        sorted_indices = sortperm(performance[end])

        threshold[i] = sorted_performance[number_of_seeds]
        pf[i] = if threshold[i] <= 0
            mean(performance[end] .<= 0)
        else
            mean(performance[end] .<= threshold[i])
        end

        ## Std MC coefficient of variation
        if i == 1
            cov[i] = sqrt((pf[i] - pf[i]^2) / sim.n) / pf[i]
            ## MarkovChain coefficient of variation
        else
            Iᵢ = reshape(
                performance[end] .< max(threshold[i], 0), samples_per_seed, number_of_seeds
            )
            cov[i] = estimate_cov(Iᵢ, pf[i], sim.n)
        end

        ## Break the loop
        if threshold[i] <= 0 || i == sim.levels
            break
        end

        nextsamples, nextperformance = nextlevelsamples(
            samples[end][sorted_indices[1:number_of_seeds], :],
            sorted_performance[1:number_of_seeds],
            threshold[i],
            models,
            performancefunction,
            inputs,
            sim,
        )

        push!(samples, nextsamples)
        push!(performance, nextperformance)
    end

    # add level to each dataframe
    for i in eachindex(samples)
        samples[i][!, :level] .= i
    end
    # merge (vcat) all samples
    samples = reduce(vcat, samples)

    pf = prod(pf)
    cov = sqrt(sum((cov .^ 2)))

    return pf, cov, samples
end

function nextlevelsamples(
    samples::DataFrame,
    performance::Vector{<:Real},
    threshold::Real,
    models::Union{Array{<:UQModel},UQModel},
    performancefunction::Function,
    inputs::Union{Array{T},T} where {T<:UQInput},
    sim::SubSetSimulation,
)
    nextlevelsamples = [samples]
    nextlevelperformance = [performance]

    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    rvs = names(random_inputs)

    number_of_chains = nrow(samples)
    samples_per_chain = Int64(floor(sim.n / number_of_chains))

    d = length(rvs)
    Φ = MvNormal(Diagonal(Matrix{Float64}(I, d, d)))

    # Modified metropolis hastings to generate samples for the next intermediate failure region
    for _ in 1:samples_per_chain
        chainsamples = copy(nextlevelsamples[end])
        chainperformance = copy(nextlevelperformance[end])

        to_standard_normal_space!(inputs, chainsamples)

        θ = Matrix{Float64}(chainsamples[:, rvs])

        ξ = θ + rand(sim.proposal, size(θ)...)
        α = pdf(Φ, transpose(ξ)) ./ pdf(Φ, transpose(θ))

        α_accept = α .>= rand(size(α)...)
        chainsamples[α_accept, rvs] = ξ[α_accept, :]

        to_physical_space!(inputs, chainsamples)

        ## Evaluating model just for new samples
        α_accept_indices = findall(α_accept)
        if any(α_accept)
            new_samples = chainsamples[α_accept, :]
            evaluate!(models, new_samples)

            new_samplesperformance = performancefunction(new_samples)
            reject = new_samplesperformance .> threshold

            new_samples[reject, :] = nextlevelsamples[end][α_accept_indices[reject], :]
            new_samplesperformance[reject] = nextlevelperformance[end][α_accept_indices[reject]]

            chainsamples[α_accept_indices, :] = new_samples
            chainperformance[α_accept_indices] = new_samplesperformance
        end

        # chainsamples = chainsamples[α_accept_indices, :]
        # chainperformance = chainperformance[α_accept_indices]

        push!(nextlevelsamples, chainsamples)
        push!(nextlevelperformance, chainperformance)
    end

    # reduce and discard seeds
    nextlevelsamples = reduce(vcat, nextlevelsamples[2:end])
    nextlevelperformance = reduce(vcat, nextlevelperformance[2:end])

    return nextlevelsamples, nextlevelperformance
end

function nextlevelsamples(
    samples::DataFrame,
    performance::Vector{<:Real},
    threshold::Real,
    models::Union{Array{<:UQModel},UQModel},
    performancefunction::Function,
    inputs::Union{Array{T},T} where {T<:UQInput},
    sim::SubSetInfinity,
)
    samples_per_seed = Int64(floor(sim.n / nrow(samples)))

    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    rvs = names(random_inputs)

    to_standard_normal_space!(inputs, samples)

    samples = repeat(samples, samples_per_seed)
    performance = repeat(performance, samples_per_seed)

    means = Matrix{Float64}(samples[:, rvs]) .* sqrt(1 - sim.s^2)

    nextlevelsamples = copy(samples)
    nextlevelsamples[:, rvs] = randn(size(means)) .* sim.s .+ means

    to_physical_space!(inputs, nextlevelsamples)

    evaluate!(models, nextlevelsamples)

    nextlevelperformance = performancefunction(nextlevelsamples)

    reject = nextlevelperformance .> threshold

    nextlevelsamples[reject, :] = samples[reject, :]
    nextlevelperformance[reject] = performance[reject]

    return nextlevelsamples, nextlevelperformance
end

"""
	estimate_cov(Iᵢ::AbstractMatrix, pf::Float64, n::Int64)

Evaluates the coefficient of variation at a subset simulation level.

Reference: Au & Beck, (2001), 'Estimation of small failure probabilities in high dimensions by subset simulation'
"""
function estimate_cov(Iᵢ::AbstractMatrix, pf::Float64, n::Int64)
    Ns, Nc = size(Iᵢ) # Number of samples per seed, number of seeds
    rᵢ = zeros(Ns - 1)
    # Eq 29 - covariance vector between indicator(l) and indicator(l+k) -> ri
    for k in 1:(Ns - 1)
        for j in 1:Nc, l in 1:(Ns - k)
            rᵢ[k] += Iᵢ[l, j] * Iᵢ[l + k, j]
        end
        rᵢ[k] /= (n - k * Nc) - pf^2
    end
    # Eq 25 - correlation coefficient vector ρ
    ρ = rᵢ / pf * (1 - pf)
    # Eq 27 - γᵢ Bernoulli coefficient
    γᵢ = 2 * sum([(1 - k * Nc / n) * ρ[k] for k in 1:(Ns - 1)])
    # Eq 28 - i-level coefficient of variation
    δᵢ = sqrt((1 - pf) / (pf * n) * (1 + γᵢ))
    return δᵢ
end
