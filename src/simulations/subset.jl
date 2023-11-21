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
        if skewness(proposal) != 0.0
            error("proposal must be a symmetric distribution")
        elseif mean(proposal) != median(proposal)
            error("proposal must be a symmetric distribution")
        elseif mean(proposal) != 0
            error("proposal must be centered in 0")
        elseif var(proposal) ≥ 2
            @warn "A proposal pdf with large variance (≥ 2) can be inefficient."
        end
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
        (0 <= s <= 1) || error("standard deviation must be between 0.0 and 1.0")
        return new(n, target, levels, s)
    end
end

function sample(inputs::Vector{<:UQInput}, sim::AbstractSubSetSimulation)
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

"""
    SubSetInfinityAdaptive(n::Integer, target::Float64, levels::Integer, Na::Integer, λ::Real, s::Real)

Implementation of: Papaioannou, Iason, et al. "MCMC algorithms for subset simulation." Probabilistic Engineering Mechanics 41 (2015): 89-103

Defines the properties of a Subset-∞ adaptive where `n` is the number of initial samples,
`target` is the target probability of failure at each level, `levels` is the maximum number
of levels and `λ` (λ = 1 recommended) is the initial scaling parameter and `Na` is the number of
times to update `λ` per subset level (number of partitions of seeds). The initial variance of the proposal distribution is `λ`.


Idea behind this algorithm is to adaptively select the correlation parameter of `s`
at each intermediate level, by simulating a subset N_a of the chains
(which must be choosen without replacement at random) and modifying the acceptance rate towards the optiming
α_star = 0.44

# Constructors
* `SubSetInfinityAdaptive(n::Integer, target::Float64, levels::Integer, Na::Integer)`   (default: λ = s = 1)
* `SubSetInfinityAdaptive(n::Integer, target::Float64, levels::Integer, Na::Integer, λ::Real)` (λ = s)
* `SubSetInfinityAdaptive(n::Integer, target::Float64, levels::Integer, Na::Integer, λ::Real, s::Real)`

# Examples

```jldoctest
julia> SubSetInfinityAdaptive(200, 0.1, 10, 2)
SubSetInfinityAdaptive(200, 0.1, 10, 2, 1, 1)
```

# References

[papaioannou2015mcmc](@cite)

[chan2022adaptive](@cite)
"""
mutable struct SubSetInfinityAdaptive <: AbstractSubSetSimulation
    n::Integer
    target::Float64
    levels::Integer
    Na::Integer
    λ::Real
    s::Real

    function SubSetInfinityAdaptive(
        n::Integer, target::Float64, levels::Integer, Na::Integer, λ::Real, s::Real
    )
        number_of_seeds = Int64(max(1, ceil(n * target)))

        (Na <= number_of_seeds) ||
            error("Number of partitions Na must be less than `n` * `target`")
        (mod(number_of_seeds, Na) == 0) ||
            error("Number of partitions Na must be a multiple of `n` * `target`")
        (0 <= λ <= 1) || error(
            "Scaling parameter must be between 0.0 and 1.0. A good initial choice is 1.0",
        )
        (0 <= s <= 1) || error("standard deviation must be between 0.0 and 1.0")
        return new(n, target, levels, Na, λ, s)
    end
end

function SubSetInfinityAdaptive(
    n::Integer, target::Float64, levels::Integer, Na::Integer, λ::Real
)
    return SubSetInfinityAdaptive(n, target, levels, Na, λ, λ)
end
function SubSetInfinityAdaptive(n::Integer, target::Float64, levels::Integer, Na::Integer)
    return SubSetInfinityAdaptive(n, target, levels, Na, 1, 1)
end

function probability_of_failure(
    models::Union{Vector{<:UQModel},UQModel},
    performancefunction::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
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
                performance[end] .<= max(threshold[i], 0), number_of_seeds, samples_per_seed
            )
            cov[i] = estimate_cov(Iᵢ, pf[i], sim.n)
        end

        @debug "Subset level $i" pf[i] threshold[i] cov[i]

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
    models::Union{Vector{<:UQModel},UQModel},
    performancefunction::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
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

    α_MCMC = zeros(samples_per_chain)
    α_ss = zeros(samples_per_chain)

    # Modified metropolis hastings to generate samples for the next intermediate failure region
    for i in 1:samples_per_chain
        chainsamples = copy(nextlevelsamples[end])
        chainperformance = copy(nextlevelperformance[end])

        to_standard_normal_space!(inputs, chainsamples)

        θ = Matrix{Float64}(chainsamples[:, rvs])

        ξ = θ + rand(sim.proposal, size(θ)...)
        α = pdf(Φ, transpose(ξ)) ./ pdf(Φ, transpose(θ))

        α_accept = α .>= rand(size(α)...)
        chainsamples[α_accept, rvs] = ξ[α_accept, :]

        α_MCMC[i] = mean(α_accept)

        to_physical_space!(inputs, chainsamples)

        ## Evaluating model just for new samples
        α_accept_indices = findall(α_accept)
        if any(α_accept)
            new_samples = chainsamples[α_accept, :]
            evaluate!(models, new_samples)

            new_samplesperformance = performancefunction(new_samples)
            reject = new_samplesperformance .> threshold

            α_ss[i] = 1 - mean(reject)

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

    @debug "acceptance rate MCMC" mean(α_MCMC)
    @debug "acceptance rate subset" mean(α_ss)

    # reduce and discard seeds
    nextlevelsamples = reduce(vcat, nextlevelsamples[2:end])
    nextlevelperformance = reduce(vcat, nextlevelperformance[2:end])

    return nextlevelsamples, nextlevelperformance
end

function nextlevelsamples(
    samples::DataFrame,
    performance::Vector{<:Real},
    threshold::Real,
    models::Union{Vector{<:UQModel},UQModel},
    performancefunction::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
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

    α = 1 - mean(reject)
    @debug "acceptance rate" α

    nextlevelsamples[reject, :] = samples[reject, :]
    nextlevelperformance[reject] = performance[reject]

    return nextlevelsamples, nextlevelperformance
end

function nextlevelsamples(
    samples::DataFrame,
    performance::Vector{<:Real},
    threshold::Real,
    models::Union{Vector{<:UQModel},UQModel},
    performancefunction::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
    sim::SubSetInfinityAdaptive,
)
    next_samples = DataFrame[]
    next_performance = Vector{eltype(performance)}[]

    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    rvs = names(random_inputs)

    a_star = 0.44 # Optimal acceptance rate, so say Papaioannou, I., et. al.

    Ns = length(performance) # Number of seeds
    permutation = shuffle(1:Ns)

    samples[:, :] = samples[permutation, :]
    performance[:] = performance[permutation]

    samples_per_seed = Integer(sim.n / Ns)
    number_of_batches = Integer(Ns / sim.Na)

    λ = sim.λ

    α = 0.0
    for i in 1:number_of_batches
        σ = min(1, λ * sim.s)
        ρ = sqrt(1 - σ^2)

        batch_samples = Vector{DataFrame}(undef, samples_per_seed)
        batch_performance = Vector{Vector{eltype(performance)}}(undef, samples_per_seed)

        batch_samples[1] = samples[((i - 1) * sim.Na + 1):(i * sim.Na), :]
        batch_performance[1] = performance[((i - 1) * sim.Na + 1):(i * sim.Na)]

        acceptance_rate = 0.0

        for k in 2:samples_per_seed
            batch_samples[k] = copy(batch_samples[k - 1])

            # to_standard_normal_space!(inputs, batch_samples[k])

            μ = Matrix{Float64}(batch_samples[k][:, rvs]) .* ρ
            batch_samples[k][:, rvs] = randn(size(μ)) .* σ .+ μ

            # to_physical_space!(inputs, batch_samples[k])

            evaluate!(models, batch_samples[k])

            batch_performance[k] = performancefunction(batch_samples[k])

            rejected = batch_performance[k] .> threshold
            batch_samples[k][rejected, rvs] = batch_samples[k - 1][rejected, rvs]
            batch_performance[k][rejected] = batch_performance[k - 1][rejected]

            acceptance_rate += (1 - mean(rejected))
        end

        ζ = 1 / sqrt(i)
        λ = λ * exp(ζ * (acceptance_rate / sim.Na - a_star))

        push!(next_samples, reduce(vcat, batch_samples))
        push!(next_performance, reduce(vcat, batch_performance))
    end

    @info "aCs" λ α / number_of_batches

    sim.λ = λ

    next_samples = reduce(vcat, next_samples)
    next_performance = reduce(vcat, next_performance)

    return next_samples, next_performance
end

"""
	estimate_cov(Iᵢ::AbstractMatrix, pf::Float64, n::Int64)

Evaluates the coefficient of variation at a subset simulation level.

Reference: Au & Beck, (2001), 'Estimation of small failure probabilities in high dimensions by subset simulation'
"""
function estimate_cov(Iᵢ::AbstractMatrix, pf::Float64, n::Int64)
    Nc, Ns = size(Iᵢ) # Number of samples per seed, number of seeds
    rᵢ = zeros(Ns - 1)
    # Eq 29 - covariance vector between indicator(l) and indicator(l+k) -> ri
    for k in 1:(Ns - 1)
        for l in 1:(Ns - k)
            rᵢ[k] += dot(Iᵢ[:, l], Iᵢ[:, l + k])
        end
        rᵢ[k] = rᵢ[k] * 1 / (n - k * Nc) - pf^2
    end
    # Eq 25 - correlation coefficient vector ρ
    ρ = rᵢ / (pf * (1 - pf))
    # Eq 27 - γᵢ Bernoulli coefficient
    γᵢ = 2 * sum([(1 - k * Nc / n) * ρ[k] for k in 1:(Ns - 1)])
    # Eq 28 - i-level coefficient of variation
    δᵢ = sqrt((1 - pf) / (pf * n) * (1 + γᵢ))
    return δᵢ
end
