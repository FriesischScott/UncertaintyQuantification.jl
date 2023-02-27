mutable struct LineSampling
    lines::Integer
    points::Vector{<:Real}
    direction::NamedTuple

    function LineSampling(
        lines::Integer,
        points::Vector{<:Real}=collect(0.5:0.5:5),
        direction::NamedTuple=NamedTuple(),
    )
        return new(lines, points, direction)
    end
end

function sample(inputs::Vector{<:UQInput}, sim::LineSampling)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)

    n_rv = count_rvs(random_inputs)
    n_samples = length(sim.points) * sim.lines
    rv_names = names(random_inputs)

    α = map(n -> sim.direction[n], rv_names)
    α /= norm(α)

    θ = rand(Normal(), n_rv, sim.lines)

    θ = θ - α * (α' * θ)
    θ = repeat(θ; outer=[length(sim.points), 1])

    θ = θ[:] + repeat(α * sim.points'; outer=[1, sim.lines])[:]

    samples = transpose(reshape(θ, n_rv, n_samples))
    samples = DataFrame(rv_names .=> eachcol(samples))

    if !isempty(deterministic_inputs)
        samples = hcat(samples, sample(deterministic_inputs, n_samples))
    end

    to_physical_space!(inputs, samples)

    return samples
end
