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

#= function probability_of_failure(
    models::Union{Array{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Array{<:UQInput},UQInput},
    sim::SubSetSimulation,
)

    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)

    samples = sample(inputs, MonteCarlo(sim.n))

    for m in models
        evaluate!(m, samples)
    end

    g_subset = performance(samples)

    number_of_chains = maximum(1, sim.n * sim.target)
    samples_per_chain = floor(sim.n / number_of_chains)

    threshold = zeros(sim.levels, 1)
    pf = zeros(sim.levels, 1)
    rejection_rates = zeros(sim.levels, 1)
    cov = zeros(sim.level, 1)

    for i ∈ 1:sim.levels
        sorted_performance = sort(g_subset)
        sorted_indices = sortperm(g_subset)

        threshold[i] = sorted_performance[number_of_chains]

        pf[i] = sum(g_subset .<= threshold[i] <= 0 ? 0 : threshold[i]) / sim.n
        cov[i] = subset_cov(sim, i, sorted_performance, threshold[i], pf[i])

        @printf("[SubSet] Estimated pf: %f", prod(pf[1:i]))

        (threshold[i] <= 0) && break

        seeds = samples[sorted_indices[1:number_of_chains], :]

    end

end =#

#= function subset_cov(sim::SubSetSimulation, level, g, threshold, pf)
    if level == 1
        return sqrt((1 - pf) / (pf * sim.n))
    end

    nc = maximum(1, sim.n * sim.target)
    sc = floor(sim.n / number_of_chains)

    ind = reshape(g[end - nc * sc + 1:end], [], sc) .< threshold

    correlation = zeros(nc, sc)

    for i ∈ 1:nc
        v = ind[i, :]
        for δk ∈ 0:sc - 1
            v1 = v[1:end - δk]
            v2 = v[1 + δk:end]
            correlation[δk + 1, i] = (1 / length(v1)) * sum(v1 .* v2)
        end
    end

    corr = sum(correlation, 2) / nc - pf^2

    ρ = corr ./ corr[1]
    γ = 2 * sum((1 - (1:sc - 1) * nc / sim.n) .* transpose(ρ[1:sc - 1]))

    sqrt((1 - pf) / (pf * sim.n) * (1 + γ))
end =#