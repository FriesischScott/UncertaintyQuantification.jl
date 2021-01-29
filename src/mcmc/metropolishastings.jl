struct MarkovChain
    targets
    proposals
    samples::Array{DataFrame,1}
end

function assemblechains(inputs, ξ, seeds)
    proposals = []
    targets = []

    # parse joint distributions
    for jd ∈ filter(i -> isa(i, JointDistribution), inputs)
        dim = dimensions(jd)
        a = ones(dim, 1) * -ξ
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
        a = ones(n, 1) * -ξ
        b = -1 .* a

        marginals = RandomVariable.(Uniform.(a, b), names(rvs))
        copula = GaussianCopula(Matrix{Float64}(I(n)))

        push!(proposals, JointDistribution(marginals, copula))
        push!(targets, JointDistribution(rvs, copula))
    end

    return MarkovChain(targets, proposals, [seeds])
end

function metropolishastings(mc::MarkovChain, n::Integer=1)

    m = size(mc.samples[1], 1)

    for i ∈ 1:n
        chain_samples = DataFrame()
        for j ∈ eachindex(mc.proposals)
            p = mc.proposals[j]
            t = mc.targets[j]

            dim = dimensions(p)
            Φ = MvNormal(dim, 1.0)

            U = sample(p, m)
            to_standard_normal_space!(p, U)
            U = convert(Matrix, U)

            U_last = mc.samples[end][:, names(t)]
            to_standard_normal_space!(t, U_last)
            U_last = convert(Matrix, U_last)

            pdf_i = prod(pdf(Φ, (U .+ U_last)'), dims=2)
            pdf_0 = prod(pdf(Φ, U_last'), dims=2)

            α = pdf_i ./ pdf_0
            accept = α .>= rand(size(α)...)

            samples = DataFrame(U_last + U .* accept)
            rename!(samples, names(p))

            to_physical_space!(t, samples)

            chain_samples = hcat(chain_samples, samples)
        end

        push!(mc.samples, chain_samples)
    end

    return mc
end