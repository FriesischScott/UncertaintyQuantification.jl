mutable struct LineSampling
    lines::Int64
    points::Vector{Float64}
    direction::DataFrame

    function LineSampling(
        lines::Int64,
        points::Vector{Float64} = collect(0.5:0.5:5),
        direction::DataFrame = DataFrame()
    )
        new(lines, points, direction)
    end
end

function sample(inputs::Array{<:UQInput}, sim::LineSampling)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)

    n_rv = 0

    for i in random_inputs
        if isa(i, RandomVariable)
            n_rv += 1
        elseif isa(i, RandomVariableSet)
            n_rv += length(i.members)
        end
    end

    α = vec(Matrix(sim.direction))
    α /= norm(α)

    θ = rand(Normal(), n_rv, sim.lines)

    θ = θ - α * (α' * θ)
    θ = repeat(θ, outer=[length(sim.points), 1])

    θ = θ[:] + repeat(α * sim.points', outer=[1, sim.lines])[:]

    samples = DataFrame(reshape(θ, n_rv, length(sim.points) * sim.lines)')

    random_names = names(random_inputs)

    names!(samples, random_names)

    if !isempty(deterministic_inputs)
        samples = hcat(samples, sample(deterministic_inputs, length(sim.points) * sim.lines))
    end

    to_physical_space!(inputs, samples)

    return samples
end
