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

    # Evaluate performance function
    p = reshape(performance(samples), length(sim.points), sim.lines)'

    # Find roots using spline interpolation for each line
    β = [rootinterpolation(sim.points, p[i, :], i) for i in 1:(sim.lines)]

    # Determine pf along lines
    ξ = cdf.(Normal(), -β)

    # Estimators for pf and variance
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
    normalize!(α)

    # Start a line from origin parallel to 𝜶, determine distance 𝛽
    θ₀ = sim.points * α'
    samples = DataFrame(rv_names .=> eachcol(θ₀))
    evaluate!(models, samples)

    β⁺ = rootinterpolation(sim.points, performance(samples))
    βᵢ = copy(β⁺)

    if isinf(β⁺)
        error("No root found on initial line")
    end

    # Generate samples in standard normal space
    θ = rand(Normal(), sim.lines, n_rv)

    # Project samples onto hyperplane orthogonal to 𝜶
    θₚ = θ - (θ * α) * α'

    # Find the sample with smallest norm
    idx = argmin(norm.(eachrow(θ)))

    # Keep track of processed indices
    notprocessed = collect(1:(sim.lines))

    # Vector of β
    β = zeros(sim.lines)

    # Loop over lines
    for i in 1:(sim.lines)
        # Calculate distance
        θᵢ = θₚ[idx, :]

        # Limit-state function along line
        f = β -> begin
            sample = DataFrame(rv_names .=> eachrow(θᵢ .+ α * β))
            evaluate!(models, sample)
            append!(samples, sample)
            return performance(sample)
        end

        β[i] = newtonraphson(βᵢ, f, sim.stepsize, sim.tolerance, sim.maxiterations)

        # Update starting point for next iteration
        if isinf(β[i])
            # Find root using interpolation if iteration does not converge
            sample = DataFrame(rv_names .=> eachcol((θᵢ .+ α * sim.points')'))
            evaluate!(models, sample)
            append!(samples, sample)
            β[i] = rootinterpolation(sim.points, performance(sample),i)
        end

        if isfinite(β[i]) βᵢ = β[i] end

        # Check if distance is smaller than previous distance
        if β[i] + 1e-6 < β⁺
            # Update β
            β⁺ = βᵢ

            # Update α
            α = normalize(θᵢ + α * βᵢ)

            # Project remaining samples onto new base
            θₚ[notprocessed, :] = θ[notprocessed, :] - (θ[notprocessed, :] * α) * α'
        end

        # Remove processed line from list
        filter!(x -> x != idx, notprocessed)

        if i !== sim.lines
            # Find next line
            idx = argmin(norm.(eachrow(θₚ[notprocessed, :] .- θᵢ')))
            idx = notprocessed[idx]
        end
    end

    pf = mean(cdf.(Normal(), -β))
    variance = var(cdf.(Normal(), -β)) / sim.lines

    to_physical_space!(inputs, samples)

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
