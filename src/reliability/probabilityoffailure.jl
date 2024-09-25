function probability_of_failure(
    models::Union{Vector{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
    sim::AbstractMonteCarlo,
)
    if isimprecise(inputs)
        error("You must use DoubleLoop or RandomSlicing with imprecise inputs.")
    end

    samples = sample(inputs, sim)
    evaluate!(models, samples)

    # Probability of failure
    pf = sum(performance(samples) .< 0) / sim.n

    variance = (pf - pf^2) / sim.n

    return pf, sqrt(variance), samples
end

function probability_of_failure(
    models::Union{Vector{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
    sim::LineSampling,
)
    if isimprecise(inputs)
        error("You must use DoubleLoop or RandomSlicing with imprecise inputs.")
    end

    if isempty(sim.direction)
        sim.direction = gradient_in_standard_normal_space(
            [models..., Model(x -> -1 * performance(x), :performance)],
            inputs,
            sns_zero_point(inputs),
            :performance,
        )
    end

    samples = sample(inputs, sim)
    evaluate!(models, samples)

    p = reshape(performance(samples), length(sim.points), sim.lines)

    ϕ = Normal()
    ξ = zeros(sim.lines)
    x = median(sim.points)
    for i in 1:(sim.lines)
        if all(p[:, i] .< 0)
            ξ[i] = 1.0
            @warn "All samples for line $i are inside the failure domain"
            continue
        elseif all(p[:, i] .> 0)
            ξ[i] = 0.0
            @warn "All samples for line $i are outside the failure domain"
            continue
        end
        spl = Spline1D(sim.points, p[:, i])
        try
            root = Dierckx.roots(spl)[1]
            ξ[i] = cdf.(ϕ, -root)
        catch e
            @warn "Intersection with failure domain not found for line $i ($e)"
        end
    end

    pf = mean(ξ)
    variance = var(ξ) / sim.lines

    return pf, sqrt(variance), samples
end

function probability_of_failure(
    models::Union{Vector{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
    sim::AdvancedLineSampling,
)
    if isimprecise(inputs)
        error("You must use DoubleLoop or RandomSlicing with imprecise inputs.")
    end

    if isempty(sim.direction)
        sim.direction = gradient_in_standard_normal_space(
            [models..., Model(x -> -1 * performance(x), :performance)],
            inputs,
            sns_zero_point(inputs),
            :performance,
        )
    end

    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)

    n_rv = count_rvs(random_inputs)
    rv_names = names(random_inputs)

    # Get the important direction 𝜶
    α = map(n -> sim.direction[n], rv_names)
    α /= norm(α)

    # Start a line from origin parallel to 𝜶, determine distance 𝛽
    θ₀ = α * sim.points'
    samples = DataFrame(rv_names .=> eachcol(θ₀'))
    β⁺ = splinefit(performance, samples, sim)
    βᵢ = copy(β⁺)

    if isinf(β⁺)
        @warn "No root found on initial line"
        return nothing
    end

    # Generate samples in standard normal space
    θ = rand(Normal(), n_rv, sim.lines)

    # Project samples onto hyperplane orthogonal to 𝜶
    θₚ = θ - α * (α' * θ)

    # Find the sample with smallest norm
    idx = argmin(norm.(eachcol(θ)))

    # Keep track of processed indices
    notprocessed = collect(1:sim.lines)

    # Vector of β
    β = zeros(sim.lines)

    # Loop over lines
    for i in 1:sim.lines
        # Calculate distance
        θᵢ = θₚ[:,idx]

        # Limit-state function along line
        f = β -> performance(DataFrame(rv_names .=> eachcol((θᵢ .+ α * β)')))
        β[i], x = newtonraphson(βᵢ, f, sim)

        append!(samples, DataFrame(rv_names .=> eachcol((θᵢ .+ α * x')')))

        # Update starting point for next iteration
        if isfinite(β[i]) βᵢ = β[i] end

        # Check if distance is smaller than previous distance
        if β[i]+1e-6 < β⁺
            # Update β
            β⁺ = βᵢ

            # Update α
            α = normalize(θᵢ + α*βᵢ)

            # Project remaining samples onto new base
            θₚ[:,notprocessed] = θ[:,notprocessed] - α * (α' * θ[:,notprocessed])
        end

        # Remove processed line from list
        filter!(x -> x != idx, notprocessed)

        if i !== sim.lines
            # Find next line
            idx = argmin(norm.(eachcol(θₚ[:,notprocessed]  .- θᵢ)))
            idx = notprocessed[idx]
        end
    end

    pf = mean(cdf.(Normal(), -β))
    variance = var(cdf.(Normal(), -β)) / sim.lines

    return pf, sqrt(variance), samples
end

function probability_of_failure(
    models::Union{Vector{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
    sim::ImportanceSampling,
)
    if isimprecise(inputs)
        error("You must use DoubleLoop or RandomSlicing with imprecise inputs.")
    end

    if isempty(sim.dp) || isempty(sim.α)
        _, β, dp, α = probability_of_failure(models, performance, inputs, FORM())
        sim.β = β
        sim.dp = dp
        sim.α = α
    end

    samples, weights = sample(inputs, sim)
    evaluate!(models, samples)

    # Probability of failure
    weighted_failures = (performance(samples) .< 0) .* weights
    pf = sum(weighted_failures) / sim.n

    variance = ((sum(weighted_failures .* weights) / sim.n) - pf^2) / sim.n

    return pf, sqrt(variance), samples
end

# Allow to calculate the pf using only a performance function but no model
function probability_of_failure(
    performance::Function, inputs::Union{Vector{<:UQInput},UQInput}, sim::Any
)
    return probability_of_failure(UQModel[], performance, wrap(inputs), sim)
end
