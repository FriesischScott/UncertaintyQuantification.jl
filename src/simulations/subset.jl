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

function probability_of_failure(models::Union{Array{<:UQModel},UQModel},performancefunction::Function,inputs::Union{Array{T},T} where {T<:UQInput},sim::SubSetSimulation)

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

        ## Std MC covariance
        if i == 1
            cov[i] = sqrt((pf[i] - pf[i]^2) / sim.n) / pf[i]
        end

        nextlevelsamples = [samples[end][sorted_indices[1:number_of_chains], :]]
        nextlevelperformance = [sorted_performance[1:number_of_chains]]

        # Modified metropolis hastings to generate samples for the next intermediate failure region
        for c in 1:samples_per_chain
            chainsamples = copy(nextlevelsamples[end])

            to_standard_normal_space!(inputs, chainsamples)

            chainsamples[:, rvs], α_accept = candidatesamples(Matrix{Float64}(chainsamples[:, rvs]), sim.proposal)

            to_physical_space!(inputs, chainsamples)

            ## Evaluating model just for new samples
            α_accept_indices = findall(x -> x == true, α_accept)
            if length(α_accept_indices) != 0
                to_eval = chainsamples[α_accept, :]
                evaluate!(models, to_eval)

                to_evalperformance = performancefunction(to_eval)
                performance_accept = to_evalperformance .< threshold[i]
                chainsamples = copy(nextlevelsamples[end])
                chainperformance = copy(nextlevelperformance[end])
                chainsamples[α_accept_indices[performance_accept], :] = to_eval[performance_accept, :]
                chainperformance[α_accept_indices[performance_accept]] = to_evalperformance[performance_accept]
            else
                chainsamples = copy(nextlevelsamples[end])
                chainperformance = copy(nextlevelperformance[end])
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

        if i > 1
            cov[i] = chaincovi(nextlevelperformance, number_of_chains, samples_per_chain, threshold[i], pf[i], sim.n)
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
    # merge (vcat) all saÊmples
    samples = reduce(vcat, samples)

    pf = prod(pf)
    cov_pf = sqrt(sum((cov .^ 2)))

    return pf, samples, cov_pf
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

function chaincovi(nextlevelperformance::Vector{Float64}, number_of_chains::Int64, samples_per_chain::Int64, threshold::Float64, pf::Float64, sim_n::Int64)
    # Building indicator Matrix (Boolean)
    indicator = reshape(nextlevelperformance, number_of_chains, samples_per_chain)
    indicator = transpose(indicator .< max(threshold, 0))
    ```Estimation of small failure probabilities in high dimensions by subset simulation
    ```
    # Eq 29 - covariance vecotr between indicator(l) and indicator(l+k) -> ri
    ri = zeros(1, samples_per_chain)
    for k in 1:samples_per_chain
        for j in 1:number_of_chains
            for l in 1:(Int(samples_per_chain - (k - 1)))
                ri[k] = ri[k] + indicator[l, j] * indicator[l + k - 1, j]
            end
        end
        ri[k] = ri[k] / (sim_n - (k - 1) * number_of_chains) - pf^2
    end
    # Eq 25 - correlation coefficient vector ρ
    ρ = ri / ri[1]
    # Eq 27 - γ_i Bernoulli coefficient 
    γ_i = 0
    for k in 1:(samples_per_chain - 1)
        γ_i = γ_i + (1 - k * number_of_chains / sim_n) * ρ[k]
    end
    γ_i = 2 * γ_i
    #Eq 28 - i-level coefficient of variation (Metropolis Markov Chain)
    cov_i = sqrt((1 - pf) / (pf * sim_n) * (1 + γ_i))
    return cov_i
end

