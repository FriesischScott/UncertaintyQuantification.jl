function sobolindices(
    models::Array{<:UQModel},
    inputs::Array{<:UQInput},
    output::Symbol,
    sim::MonteCarlo
)

    A = sample(inputs, sim.n)
    B = sample(inputs, sim.n)

    random_names = filter(i -> isa(i, RandomUQInput), inputs) |> names

    for m ∈ models
        evaluate!(m, A) 
        evaluate!(m, B)
    end

    fA = A[:, output]
    fB = B[:, output]
    f0 = 1/sim.n * sum(fA .* fB) # given estimator for (f_0)^2

    VY = 1/sim.n * sum(fA.^2) - f0
    Vi = zeros(length(random_names))

    for (i, name) ∈ enumerate(random_names)
        ABi = copy(A)
        ABi[:, name] = B[:, name]

        for m ∈ models
            evaluate!(m, ABi)
        end

        Vi[i] = 1/(2*sim.n) * sum((fA .- ABi[:, output]).^2)
    end

    (;zip(random_names, Vi ./ VY)...)
end
