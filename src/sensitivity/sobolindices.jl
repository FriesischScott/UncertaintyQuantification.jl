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

    indices = Dict()

    A = samples[1:(sim.n), :]
    B = samples[(sim.n + 1):end, :]

    fA = Matrix(A[!, outputs])
    fB = Matrix(B[!, outputs])
    fA = fA .- mean(eachrow(fA))'
    fB = fB .- mean(eachrow(fB))'
    VY = var.(collect(eachcol(vcat(fA, fB))))
    Si = zeros(length(outputs), length(random_names), 2)
    STi = zeros(length(outputs), length(random_names), 2)

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
            Si[j, i, 1] = first_order(ABi[:, qty])
            bs = bootstrap(first_order, ABi[:, qty], BasicSampling(1000))
            Si[j, i, 2] = stderror(bs)[1]
            # Total effects
            STi[j, i, 1] = total_effect(ABi[:, qty])
            bs = bootstrap(total_effect, ABi[:, qty], BasicSampling(1000))
            STi[j, i, 2] = stderror(bs)[1]
        end
    end

    indices = Dict()
    for (i, output) in enumerate(outputs)
        indices[output] = DataFrame(
            Variables = random_names,
            FirstOrder = Si[i, :, 1],
            FirstOrderStdError = Si[i, :, 2],
            TotalEffect = STi[i, :, 1],
            TotalEffectStdError = STi[i, :, 2]
        )
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
