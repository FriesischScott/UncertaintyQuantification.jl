function sobolindices(
    models::Array{<:UQModel},
    inputs::Array{<:UQInput},
    output::Symbol,
    sim::AbstractMonteCarloSampling
)

    sim_double_samples = @set sim.n = 2 * sim.n

    samples = sample(inputs, sim_double_samples)

    random_names = filter(i -> isa(i, RandomUQInput), inputs) |> names

    for m ∈ models
        evaluate!(m, samples)
    end
    
    A = samples[1:sim.n, :]
    B = samples[sim.n + 1:end, :]

    fA = A[:, output]
    fB = B[:, output]
    f0 = 1 / sim.n * sum(fA .* fB) # given estimator for (f_0)^2

    VY = 1 / sim.n * sum(fA.^2) - f0
    Vi = zeros(length(random_names))

    for (i, name) ∈ enumerate(random_names)
        ABi = copy(A)
        ABi[:, name] = B[:, name]

        for m ∈ models
            evaluate!(m, ABi)
        end

        Vi[i] = 1 / (2 * sim.n) * sum((fA .- ABi[:, output]).^2)
    end

    (;zip(random_names, Vi ./ VY)...)
end
