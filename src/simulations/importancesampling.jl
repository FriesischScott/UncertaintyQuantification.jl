struct ImportanceSampling
    n::Integer
    β::Real
    dp::NamedTuple
    α::NamedTuple
    function ImportanceSampling(n, β, dp, α) 
        @assert n > 0 "n must be greater than zero"
        new(n, β, dp, α)
    end
end

function sample(inputs::Vector{<:UQInput}, sim::ImportanceSampling)
    β = sim.β
    dp = DataFrame([sim.dp])
    to_standard_normal_space!(inputs, dp)
    α = -1 * collect(sim.α)

    random_names = names(filter(i -> isa(i, RandomUQInput), inputs))
    deterministic_names = names(filter(i -> isa(i, DeterministicUQInput), inputs))

    # generate [n x m] samples in SNS
    Z = randn(sim.n, length(random_names))

    # samples perpendicular to important direction
    Z -= (Z * α) * α'

    # force in direction parallel to important direction
    b = exp(-(β^2/2)) / (cdf(Normal(), -β) * sqrt(2 * π)) # mean
    v = 2(b - β) # std, c = 2
    Zforced = randn(sim.n) .* v .+ b

    Z += Zforced * α'

    weights = DataFrame(ones(sim.n, 3), [:f, :h, :w]) 
    weights.f .= pdf.(Normal(), Zforced)
    weights.h .= pdf.(Normal(b, v), Zforced)
    weights.w = weights.f ./ weights.h

    samples = DataFrame(Z, random_names)
    to_physical_space!(filter(i -> isa(i, RandomUQInput), inputs), samples)

    if !isempty(deterministic_names)
        deterministic_samples = sample(filter(i -> isa(i, DeterministicUQInput), inputs), sim.n)
        samples = hcat(deterministic_samples, samples)
    end

    return samples, weights
end