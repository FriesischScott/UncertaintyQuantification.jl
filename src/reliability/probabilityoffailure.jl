function probability_of_failure(
    models::Union{Vector{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
    sim::AbstractMonteCarlo,
)
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
    if isempty(sim.direction)
        direction = gradient_in_standard_normal_space(
            wrap(models),
            inputs,
            DataFrame(names(inputs) .=> sns_zero_point(inputs)),
            x -> -1 * performance(x),
        )
        sim_next = LineSampling(sim.lines, sim.points, direction)
    else
        sim_next = deepcopy(sim)
    end

    samples = sample(inputs, sim_next)
    evaluate!(models, samples)

    p = reshape(performance(samples), length(sim.points), sim_next.lines)

    ϕ = Normal()
    ξ = zeros(sim_next.lines)
    x = median(sim_next.points)
    for i in 1:(sim_next.lines)
        if all(p[:, i] .< 0)
            ξ[i] = 1.0
            @warn "All samples for line $i are inside the failure domain"
            continue
        elseif all(p[:, i] .> 0)
            ξ[i] = 0.0
            @warn "All samples for line $i are outside the failure domain"
            continue
        end
        spl = Spline1D(sim_next.points, p[:, i])
        try
            root = Dierckx.roots(spl)[1]
            ξ[i] = cdf.(ϕ, -root)
        catch e
            @warn "Intersection with failure domain not found for line $i ($e)"
        end
    end

    pf = mean(ξ)
    variance = var(ξ) / sim_next.lines

    return pf, sqrt(variance), samples
end

function probability_of_failure(
    models::Union{Vector{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
    sim::ImportanceSampling,
)
    if isempty(sim.dp) || isempty(sim.α)
        _, β, dp, α = probability_of_failure(models, performance, inputs, FORM())
        sim_next = ImportanceSampling(sim.n, β, dp, α; c=sim.c)
    else
        sim_next = deepcopy(sim)
    end

    samples, weights = sample(inputs, sim_next)
    evaluate!(models, samples)

    # Probability of failure
    weighted_failures = (performance(samples) .< 0) .* weights
    pf = sum(weighted_failures) / sim_next.n

    variance = ((sum(weighted_failures .* weights) / sim_next.n) - pf^2) / sim_next.n

    return pf, sqrt(variance), samples
end

# Allow to calculate the pf using only a performance function but no model
function probability_of_failure(
    performance::Function, inputs::Union{Vector{<:UQInput},UQInput}, sim::Any
)
    return probability_of_failure(UQModel[], performance, wrap(inputs), sim)
end
