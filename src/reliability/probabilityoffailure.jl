function probability_of_failure(
    models::Union{Array{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Array{<:UQInput},UQInput},
    sim::AbstractMonteCarlo,
)
    samples = sample(inputs, sim)
    evaluate!(models, samples)

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
            :performance,
        )
    end

    samples = sample(inputs, sim)
    evaluate!(models, samples)

    p = reshape(performance(samples), length(sim.points), sim.lines)

    ϕ = Normal()
    pf = 0
    roots_found = 0
    x = median(sim.points)
    for i in 1:(sim.lines)
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

# Allow to calculate the pf using only a performance function but no model
function probability_of_failure(
    performance::Function, inputs::Union{Array{<:UQInput},UQInput}, sim::Any
)
    return probability_of_failure(UQModel[], performance, inputs, sim)
end
