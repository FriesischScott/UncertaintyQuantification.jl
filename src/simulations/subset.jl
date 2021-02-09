struct SubSetSimulation
    n::Integer
    target::Float64
    levels::Integer
    parameter
end

function sample(inputs::Array{<:UQInput}, sim::SubSetSimulation)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)

    n_rv = count_rvs(random_inputs)

    samples = rand(MvNormal(n_rv, 1), sim.n) |> transpose |> DataFrame
    rename!(samples, names(random_inputs))

    to_physical_space!(random_inputs, samples)

    if !isempty(deterministic_inputs)
        samples = hcat(samples, sample(deterministic_inputs, sim.n))
    end

    return samples
end

function probability_of_failure(
    models::Union{Array{<:UQModel},UQModel},
    performancefunction::Function,
    inputs::Union{Array{T},T} where T <: UQInput,
    sim::SubSetSimulation,
)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)

    rvs = names(random_inputs)
    samples = [sample(inputs, sim)]

    for m in models
        evaluate!(m, samples[end])
    end

    performance = [performancefunction(samples[end])]

    proposals, targets = build_markov_chains(sim, inputs)

    number_of_chains = max(1, ceil(sim.n * sim.target)) |> Int64
    samples_per_chain = floor(sim.n / number_of_chains)

    threshold = zeros(sim.levels, 1)
    pf = zeros(sim.levels, 1)

    if !isempty(deterministic_inputs)
        deterministic_samples = sample(deterministic_inputs, number_of_chains)
    else
        deterministic_samples = DataFrame()
    end

    for i ∈ 1:sim.levels
        sorted_performance = sort(performance[end])
        sorted_indices = sortperm(performance[end])

        threshold[i] = sorted_performance[number_of_chains]

        pf[i] = sum(performance[end] .<= (threshold[i] <= 0 ? 0 : threshold[i])) / sim.n

        if (threshold[i] <= 0) || i == sim.levels
            break
        end

        seeds = samples[end][sorted_indices[1:number_of_chains], :]

        nextlevelsamples = [seeds]
        nextlevelperformance = [sorted_performance[1:number_of_chains]]

        for c ∈ 1:samples_per_chain
            chainsamples = metropolishastings(proposals, targets, nextlevelsamples[end])

            chainsamples = hcat(chainsamples, deterministic_samples)

            for m ∈ models
                evaluate!(m, chainsamples)
            end

            chainperformance = performancefunction(chainsamples)

            reject = chainperformance .> threshold[i]

            chainperformance[reject] = nextlevelperformance[end][reject]
            chainsamples[reject, rvs] = nextlevelsamples[end][reject, rvs]

            if c == 1
                nextlevelsamples = [chainsamples]
                nextlevelperformance = [chainperformance]
            else
                push!(nextlevelsamples, chainsamples)
                push!(nextlevelperformance, chainperformance)
            end
        end

        nextlevelsamples = reduce(vcat, nextlevelsamples)
        nextlevelperformance = reduce(vcat, nextlevelperformance)

        push!(samples, nextlevelsamples)
        push!(performance, nextlevelperformance)
    end

    # add level to each dataframe
    for i ∈ eachindex(samples)
        samples[i][!, :level] .= i
    end
    # merge (vcat) all samples
    samples = reduce(vcat, samples)

    return prod(pf[pf .> 0]), samples
end

function build_markov_chains(sim::SubSetSimulation, inputs)
    proposals = []
    targets = []

    # parse joint distributions
    for jd ∈ filter(i -> isa(i, JointDistribution), inputs)
        dim = dimensions(jd)
        a = ones(dim, 1) * -sim.parameter
        b = -1 .* a

        marginals = RandomVariable.(Uniform.(a, b), names(jd))
        copula = GaussianCopula(Matrix{Float64}(I(dim)))

        push!(targets, jd)
        push!(proposals, JointDistribution(marginals, copula))
    end

    # parse random variables
    rvs = filter(i -> isa(i, RandomVariable), inputs)
    rvs = convert(Array{RandomVariable,1}, rvs)
    n = length(rvs)

    if n > 0
        a = ones(n, 1) * -sim.parameter
        b = -1 .* a

        marginals = RandomVariable.(Uniform.(a, b), names(rvs))
        copula = GaussianCopula(Matrix{Float64}(I(n)))

        push!(proposals, JointDistribution(marginals, copula))
        push!(targets, JointDistribution(rvs, copula))
    end

    return proposals, targets
end

function metropolishastings(proposals, targets, samples::DataFrame)
    mcsamples = map((p, t) -> begin
        dim = dimensions(p)
        Φ = MvNormal(dim, 1.0)

        U = sample(p, size(samples, 1))
        to_standard_normal_space!(p, U)
        U = convert(Matrix, U)

        U_last = samples[:, names(t)]
        to_standard_normal_space!(t, U_last)
        U_last = convert(Matrix, U_last)

        pdf_i = prod(pdf(Φ, (U .+ U_last)'), dims=2)
        pdf_0 = prod(pdf(Φ, U_last'), dims=2)

        α = pdf_i ./ pdf_0
        accept = α .>= rand(size(α)...)

        s = DataFrame(names(p) .=> eachcol(U_last + U .* accept))

        to_physical_space!(t, s)

        return s
    end, proposals, targets)

    return reduce(hcat, mcsamples)
end