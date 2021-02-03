function probability_of_failure(
    models::Union{Array{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Array{<:UQInput},UQInput},
    sim::AbstractMonteCarloSampling,
)

    samples = sample(inputs, sim.n)

    # Models
    for m in models
        evaluate!(m, samples)
    end

    # Probability of failure
    pf = sum(performance(samples) .< 0) / sim.n

    return pf, samples
end

function probability_of_failure(
    models::Union{Array{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Array{<:UQInput},UQInput},
    sim::LineSampling,
)

    if isempty(sim.direction)
        sim.direction = gradient_in_standard_normal_space(
        [models..., Model(x -> -1 * performance(x), :performance)],
        inputs,
        mean(inputs),
        :performance)
    end

    samples = sample(inputs, sim)

    for m in models
        evaluate!(m, samples)
    end

    p = reshape(performance(samples), length(sim.points), sim.lines)

    ϕ = Normal()
    pf = 0
    roots_found = 0
    x = median(sim.points)
    for i = 1:sim.lines
        if all(p[:, i] .< 0)
            @warn "All samples for line $i are inside the failure domain"
            continue
        elseif all(p[:, i] .> 0)
            @warn "All samples for line $i are outside the failure domain"
            continue
        end
        spl = Spline1D(sim.points, p[:, i])
        try
            root = Dierckx.roots(spl)[1]
            pf += cdf.(ϕ, -1 * root)
            roots_found += 1
        catch e
            @warn "Intersection with failure domain not found for line $i ($e)"
        end
    end

    pf /= roots_found

    return pf, samples
end

function probability_of_failure(
    models::Union{Array{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Array{T},T} where T <: UQInput,
    sim::UncertaintyQuantification.SubSetSimulation,
)

    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)

    levelsamples = sample(inputs, sim)

    for m in models
        evaluate!(m, levelsamples)
    end

    samples = copy(levelsamples)
    samples[!, :level] .= 1

    g_subset = performance(levelsamples)

    number_of_chains = max(1, ceil(sim.n * sim.target)) |> Int64
    samples_per_chain = floor(sim.n / number_of_chains)

    threshold = zeros(sim.levels, 1)
    pf = zeros(sim.levels, 1)
    rejection_rates = zeros(sim.levels, 1)
    cov = zeros(sim.levels, 1)

    for i ∈ 1:sim.levels
        sorted_performance = sort(g_subset)
        sorted_indices = sortperm(g_subset)

        threshold[i] = sorted_performance[number_of_chains]

        pf[i] = sum(g_subset .<= (threshold[i] <= 0 ? 0 : threshold[i])) / sim.n

        if (threshold[i] <= 0) || i == sim.levels
            break
        end

        seeds = levelsamples[sorted_indices[1:number_of_chains], names(random_inputs)]

        # Markov chain / next level samples
        markovchains = assemblechains(random_inputs, sim.parameter, seeds)
        lastperformance = sorted_performance[1:number_of_chains]

        for c ∈ 1:samples_per_chain
            markovchains = metropolishastings(markovchains)

            chainsamples = copy(markovchains.samples[end])

            if length(deterministic_inputs) > 0
                chainsamples = hcat(chainsamples, sample(deterministic_inputs, number_of_chains))
            end

            for m in models
                evaluate!(m, chainsamples)
            end

            chainperformance = performance(chainsamples)

            reject = chainperformance .> threshold[i]

            chainperformance[reject] = lastperformance[reject]
            markovchains.samples[end][reject, :] = markovchains.samples[end - 1][reject, :]

            lastperformance = chainperformance

            if c == 1
                levelsamples = chainsamples
                g_subset = chainperformance
            else
                append!(levelsamples, chainsamples)
                g_subset = vcat(g_subset, chainperformance)
            end
        end
        levelsamples[!, :level] .= i + 1
        append!(samples, levelsamples)
    end

    return prod(pf[pf .> 0]), samples
end