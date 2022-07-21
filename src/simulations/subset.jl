struct SubSetSimulation
    n::Integer
    target::Float64
    levels::Integer
    proposal::Sampleable{Univariate}

    function SubSetSimulation(
        n::Integer, target::Float64, levels::Integer, proposal::Sampleable{Univariate}
    )
        skewness(proposal) != 0.0 && error("proposal must be a symmetric distribution")
        mean(proposal) != median(proposal) &&
            error("proposal must be a symmetric distribution")
        return new(n, target, levels, proposal)
    end
end

function sample(inputs::Array{<:UQInput}, sim::SubSetSimulation)
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
    sim::SubSetSimulation,
)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    rvs = names(random_inputs)

    samples = [sample(inputs, sim)]

    evaluate!(models, samples[end])

    performance = [performancefunction(samples[end])]

    number_of_chains = Int64(max(1, ceil(sim.n * sim.target)))
    samples_per_chain = Int64(floor(sim.n / number_of_chains))

    threshold = zeros(sim.levels, 1)
    pf = ones(sim.levels, 1)
    cov = zeros(sim.levels, 1)

    for i in 1:(sim.levels)
        sorted_performance = sort(performance[end])
        sorted_indices = sortperm(performance[end])

        threshold[i] = sorted_performance[number_of_chains]
        pf[i] = if threshold[i] <= 0
            mean(performance[end] .<= 0)
        else
            mean(performance[end] .<= threshold[i])
        end

        nextlevelsamples = [samples[end][sorted_indices[1:number_of_chains], :]]
        nextlevelperformance = [sorted_performance[1:number_of_chains]]

        # Modified metropolis hastings to generate samples for the next intermediate failure region
        for c in 1:samples_per_chain
            chainsamples = copy(nextlevelsamples[end])
            chainperformance = copy(nextlevelperformance[end])

            to_standard_normal_space!(inputs, chainsamples)

            chainsamples[:, rvs], α_accept = candidatesamples(Matrix{Float64}(chainsamples[:, rvs]), sim.proposal)

            to_physical_space!(inputs, chainsamples)

            ## Evaluating model just for new samples
            α_accept_indices = findall(α_accept)
            if any(α_accept)
                new_samples = chainsamples[α_accept, :]
                evaluate!(models, new_samples)

                new_samplesperformance = performancefunction(new_samples)
                performance_accept = new_samplesperformance .< threshold[i]
                chainsamples = copy(nextlevelsamples[end])
                chainperformance = copy(nextlevelperformance[end])
                chainsamples[α_accept_indices[performance_accept], :] = new_samples[performance_accept, :]
                chainperformance[α_accept_indices[performance_accept]] = new_samplesperformance[performance_accept]
            end

            if c == 1
                nextlevelsamples = [chainsamples]
                nextlevelperformance = [chainperformance]
            else
                push!(nextlevelsamples, chainsamples)
                push!(nextlevelperformance, chainperformance)
            end
        end

        nextlevelsamples = reduce(vcat, nextlevelsamples)
        nextlevelperformance = reduce(vcat, nextlevelperformance)

        push!(samples, nextlevelsamples)
        push!(performance, nextlevelperformance)

        ## Std MC coefficient of variation
        if i == 1
            cov[i] = sqrt((pf[i] - pf[i]^2) / sim.n) / pf[i]
            ## MarkovChain coefficient of variation
        else
            Iᵢ = reshape(
                nextlevelperformance .< max(threshold[i], 0),
                samples_per_chain,
                number_of_chains,
            )
            cov[i] = estimate_chain_cov(Iᵢ, pf[i], sim.n)
        end
        ## Break the loop
        if threshold[i] <= 0 || i == sim.levels
            break
        end
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

function candidatesamples(θ::AbstractMatrix, proposal::Sampleable{Univariate})
    d = size(θ, 2)
    Φ = MvNormal(Diagonal(Matrix{Float64}(I, d, d)))

    ξ = θ + rand(proposal, size(θ)...)
    α = pdf(Φ, transpose(ξ)) ./ pdf(Φ, transpose(θ))

    accept = α .>= rand(size(α)...)
    θ[accept, :] = ξ[accept, :]

    return θ, accept
end

"""
	sample(inputs::Array{<:UQInput}, n::Integer)

Evaluates coefficient of variation of each subset simulation's level.
Reference: 'Estimation of small failure probabilities in high dimensions by subset simulation' - Siu-Kui Au, James L. Beck
        - Eq 29 - covariance vecotr between indicator(l) and indicator(l+k) -> ri
        - Eq 25 - correlation coefficient vector ρ
        - Eq 27 - γᵢ Bernoulli coefficient 
        - Eq 28 - i-level coefficient of varᵢation (Metropolis Markov Chain)
"""
function estimate_chain_cov(Iᵢ, pf::Float64, n::Int64)
    (samples_per_chain, Nc) = size(Iᵢ)
    rᵢ = zeros(samples_per_chain)
    for k in 1:samples_per_chain
        for j in 1:Nc, l in 1:(samples_per_chain - (k - 1))
            rᵢ[k] = rᵢ[k] + I[l, j] * I[l + k - 1, j]
        end
        rᵢ[k] = rᵢ[k] / (n - (k - 1) * Nc) - pf^2
    end
    ρ = rᵢ / rᵢ[1]
    γᵢ = 2 * sum((1 .- ((1:(samples_per_chain - 1)) .* Nc ./ n)) .* ρ[1:(end - 1)])
    δᵢ = sqrt((1 - pf) / (pf * n) * (1 + γᵢ))
    return δᵢ
end