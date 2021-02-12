struct SubSetSimulation
    n::Integer
    target::Float64
    levels::Integer
    proposal::Sampleable{Univariate}

    function SubSetSimulation(
        n::Integer,
        target::Float64,
        levels::Integer,
        proposal::Sampleable{Univariate}
    )
        skewness(proposal) != 0.0 && error("proposal must be a symmetric distribution")
        mean(proposal) != median(proposal) && error("proposal must be a symmetric distribution")
        new(n, target, levels, proposal)
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
    inputs::Union{Array{T},T} where T <: UQInput,
    sim::SubSetSimulation,
)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    rvs = names(random_inputs)

    samples = [sample(inputs, sim)]

    evaluate!(models, samples[end])

    performance = [performancefunction(samples[end])]

    number_of_chains = max(1, ceil(sim.n * sim.target)) |> Int64
    samples_per_chain = floor(sim.n / number_of_chains)

    threshold = zeros(sim.levels, 1)
    pf = ones(sim.levels, 1)

    for i ∈ 1:sim.levels
        sorted_performance = sort(performance[end])
        sorted_indices = sortperm(performance[end])

        threshold[i] = sorted_performance[number_of_chains]
        pf[i] = threshold[i] <= 0 ? mean(performance[end] .<= 0) : mean(performance[end] .<= threshold[i])

        if threshold[i] <= 0 || i == sim.levels
            break
        end

        nextlevelsamples = [samples[end][sorted_indices[1:number_of_chains], :]]
        nextlevelperformance = [sorted_performance[1:number_of_chains]]

        # Modified metropolis hastings to generate samples for the next intermediate failure region
        for c ∈ 1:samples_per_chain
            chainsamples = copy(nextlevelsamples[end])

            to_standard_normal_space!(inputs, chainsamples)

            chainsamples[:, rvs] = candidatesamples(convert(Matrix, chainsamples[:, rvs]), sim.proposal)

            to_physical_space!(inputs, chainsamples)

            evaluate!(models, chainsamples)

            chainperformance = performancefunction(chainsamples)

            reject = chainperformance .> threshold[i]

            chainperformance[reject] = nextlevelperformance[end][reject]
            chainsamples[reject, rvs] = nextlevelsamples[end][reject, rvs]

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
    end

    # add level to each dataframe
    for i ∈ eachindex(samples)
        samples[i][!, :level] .= i
    end
    # merge (vcat) all samples
    samples = reduce(vcat, samples)

    pf = prod(pf)

    return pf, samples
end

function candidatesamples(θ::AbstractMatrix, proposal::Sampleable{Univariate})
    Φ = MvNormal(size(θ, 2), 1.0)

    ξ = θ + rand(proposal, size(θ)...)
    α = pdf(Φ, transpose(ξ)) ./ pdf(Φ, transpose(θ))

    accept = α .>= rand(size(α)...)
    θ[accept, :] = ξ[accept, :]

    return θ
end