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

    fA = zeros(sim.n, length(output))
    fB = zeros(sim.n, length(output))
    VY = Dict{Symbol,Number}()
    for (i, qty) in enumerate(output)
        fA[:, i] = A[:, qty]
        fB[:, i] = B[:, qty]

        fA[:, i] .-= mean(fA[:, i])
        fB[:, i] .-= mean(fB[:, i])

        VY[qty] = var([fA[:, i]; fB[:, i]])
    end

    Si = zeros(length(random_names), 2)
    STi = zeros(length(random_names), 2)

    for (i, name) in enumerate(random_names)
        ABi = select(A, Not(name))
        ABi[:, name] = B[:, name]

        for m in models
            evaluate!(m, ABi)
        end

        for (j, qty) in enumerate(output)
            ABi[:, qty] .-= mean(ABi[:, qty])

            first_order = x -> mean(fB[:, j] .* (x .- fA[:, j])) / VY[qty] # Saltelli 2009
            total_effect = x -> (1 / (2 * sim.n)) * sum((fA[:, j] .- x) .^ 2) / VY[qty] # Saltelli 2009

            # First order effects
            Si[i, 1] = first_order(ABi[:, qty])
            bs = bootstrap(first_order, ABi[:, qty], BasicSampling(1000))
            Si[i, 2] = stderror(bs)[1]

            # Total effects
            STi[i, 1] = total_effect(ABi[:, qty])
            bs = bootstrap(total_effect, ABi[:, qty], BasicSampling(1000))
            STi[i, 2] = stderror(bs)[1]

            indices[qty] = DataFrame()
            indices[qty].Variables = random_names
            indices[qty].FirstOrder = Si[:, 1]
            indices[qty].FirstOrderStdError = Si[:, 2]
            indices[qty].TotalEffect = STi[:, 1]
            indices[qty].TotalEffectStdError = STi[:, 2]
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
