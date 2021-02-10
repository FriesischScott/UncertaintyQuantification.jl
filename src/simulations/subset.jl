struct SubSetSimulation
    n::Integer
    target::Float64
    levels::Integer
    proposal::Sampleable{Univariate}
end

function sample(inputs::Array{<:UQInput}, sim::SubSetSimulation)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)

    n_rv = count_rvs(random_inputs)

    samples = rand(MvNormal(n_rv, 1), sim.n) |> transpose |> DataFrame
    rename!(samples, names(random_inputs))

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
    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)

    rvs = names(random_inputs)
    samples = [sample(inputs, sim)]

    for m in models
        evaluate!(m, samples[end])
    end

    performance = [performancefunction(samples[end])]

    number_of_chains = max(1, ceil(sim.n * sim.target)) |> Int64
    samples_per_chain = floor(sim.n / number_of_chains)

    threshold = zeros(sim.levels, 1)
    pf = zeros(sim.levels, 1)

    if !isempty(deterministic_inputs)
        deterministic_samples = sample(deterministic_inputs, number_of_chains)
    else
        deterministic_samples = DataFrame()
    end

    for i ∈ 1:sim.levels
        sorted_performance = sort(performance[end])
        sorted_indices = sortperm(performance[end])

        threshold[i] = sorted_performance[number_of_chains]

        pf[i] = sum(performance[end] .<= (threshold[i] <= 0 ? 0 : threshold[i])) / sim.n

        if (threshold[i] <= 0) || i == sim.levels
            break
        end

        nextlevelsamples = [samples[end][sorted_indices[1:number_of_chains], :]]
        nextlevelperformance = [sorted_performance[1:number_of_chains]]

        for c ∈ 1:samples_per_chain
            chainsamples = copy(nextlevelsamples[end])

            to_standard_normal_space!(inputs, chainsamples)

            chainsamples[:, rvs] = metropolishastings(convert(Matrix, chainsamples[:, rvs]), sim.proposal)

            to_physical_space!(inputs, chainsamples)

            for m ∈ models
                evaluate!(m, chainsamples)
            end

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

    return prod(pf[pf .> 0]), samples
end