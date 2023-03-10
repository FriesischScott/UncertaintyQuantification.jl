function sobolindices(
    models::Vector{<:UQModel},
    inputs::Vector{<:UQInput},
    outputs::Vector{Symbol},
    sim::AbstractMonteCarlo,
)
    sim_double_samples = @set sim.n = 2 * sim.n

    samples = sample(inputs, sim_double_samples)
    samples = samples[shuffle(1:size(samples, 1)), :]

    random_names = names(filter(i -> isa(i, RandomUQInput), inputs))

    evaluate!(models, samples)
    indices = Dict([
        i => DataFrame(;
            Variable=Any[],
            FirstOrder=Float64[],
            FirstOrderStdError=Float64[],
            TotalEffect=Float64[],
            TotalEffectError=Float64[],
        ) for i in outputs
    ])

    A = samples[1:(sim.n), :]
    B = samples[(sim.n + 1):end, :]

    fA = Matrix(A[:, outputs])
    fB = Matrix(B[:, outputs])
    fA = fA .- mean(fA; dims=1)
    fB = fB .- mean(fB; dims=1)

    VY = var([fA; fB]; dims=1)

    for (i, name) in enumerate(random_names)
        ABi = select(A, Not(name))
        ABi[:, name] = B[:, name]

        for m in models
            evaluate!(m, ABi)
        end

        for (j, qty) in enumerate(outputs)
            ABi[:, qty] .-= mean(ABi[:, qty])

            first_order = x -> mean(fB[:, j] .* (x .- fA[:, j])) / VY[j] # Saltelli 2009
            total_effect = x -> (1 / (2 * sim.n)) * sum((fA[:, j] .- x) .^ 2) / VY[j] # Saltelli 2009

            # First order effects
            bs = bootstrap(first_order, ABi[:, qty], BasicSampling(1000))
            # Total effects
            bs = bootstrap(total_effect, ABi[:, qty], BasicSampling(1000))
            push!(
                indices[qty],
                [name first_order(ABi[:, qty]) stderror(bs)[1] total_effect(ABi[:, qty]) stderror(
                    bs
                )[1]],
            )
        end
    end
    return indices
end

function sobolindices(pce::PolynomialChaosExpansion)
    random_names = names(pce.inputs)
    indices = DataFrame()
    indices.Variables = random_names
    indices.FirstOrder =
        [
            sum(pce.y[findall(α -> α[i] != 0 && sum(α) == α[i], pce.Ψ.α)] .^ 2) for
            i in 1:length(random_names)
        ] ./ var(pce)
    indices.TotalEffect =
        [
            sum(pce.y[findall(α -> α[i] != 0, pce.Ψ.α)] .^ 2) for
            i in 1:length(random_names)
        ] ./ var(pce)

    return indices
end
