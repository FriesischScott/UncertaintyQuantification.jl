function probability_of_failure(
    models::Union{Vector{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
    sim::AbstractMonteCarlo,
)
    @assert !isimprecise(inputs)

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
    @assert !isimprecise(inputs)

    if isempty(sim.direction)
        sim.direction = gradient_in_standard_normal_space(
            [models..., Model(x -> -1 * performance(x), :performance)],
            inputs,
            DataFrame(names(inputs) .=> mean(inputs)),
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
    sim::ImportanceSampling,
)
    @assert !isimprecise(inputs)

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
