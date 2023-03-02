function sobolindices(
    models::Vector{<:UQModel},
    inputs::Vector{<:UQInput},
    output::Vector{Symbol},
    sim::AbstractMonteCarlo,
)
    sim_double_samples = @set sim.n = 2 * sim.n

    samples = sample(inputs, sim_double_samples)
    samples = samples[shuffle(1:size(samples, 1)), :]

    random_names = names(filter(i -> isa(i, RandomUQInput), inputs))

    evaluate!(models, samples)

    indices = Dict()

    A = samples[1:(sim.n), :]
    B = samples[(sim.n + 1):end, :]

    fA = DataFrame()
    fB = DataFrame()
    VY = Dict{Symbol,Number}()
    for quantity in output
        fA[!, quantity] = A[:, quantity]
        fB[!, quantity] = B[:, quantity]

        fA[!, quantity] .-= mean(fA[!, quantity])
        fB[!, quantity] .-= mean(fB[!, quantity])

        VY[quantity] = var([fA[!, quantity]; fB[!, quantity]])
    end

    Si = zeros(length(random_names), 2)
    STi = zeros(length(random_names), 2)

    for (i, name) in enumerate(random_names)
        ABi = select(A, Not(name))
        ABi[:, name] = B[:, name]

        for m in models
            evaluate!(m, ABi)
        end

        for quantity in output
            ABi[:, quantity] .-= mean(ABi[:, quantity])

            first_order =
                x -> mean(fB[!, quantity] .* (x .- fA[!, quantity])) / VY[quantity] # Saltelli 2009
            total_effect =
                x -> (1 / (2 * sim.n)) * sum((fA[!, quantity] .- x) .^ 2) / VY[quantity] # Saltelli 2009

            # First order effects
            Si[i, 1] = first_order(ABi[:, quantity])
            bs = bootstrap(first_order, ABi[:, quantity], BasicSampling(1000))
            Si[i, 2] = stderror(bs)[1]

            # Total effects
            STi[i, 1] = total_effect(ABi[:, quantity])
            bs = bootstrap(total_effect, ABi[:, quantity], BasicSampling(1000))
            STi[i, 2] = stderror(bs)[1]

            indices[quantity] = DataFrame()
            indices[quantity].Variables = random_names
            indices[quantity].FirstOrder = Si[:, 1]
            indices[quantity].FirstOrderStdError = Si[:, 2]
            indices[quantity].TotalEffect = STi[:, 1]
            indices[quantity].TotalEffectStdError = STi[:, 2]
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
