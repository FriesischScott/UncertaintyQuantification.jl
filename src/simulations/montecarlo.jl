struct MonteCarlo <: AbstractMonteCarloSampling
    n::Integer
    MonteCarlo(n) = n > 0 ? new(n) : error("n must be greater than zero")
end

struct SobolSampling <: AbstractMonteCarloSampling
    n::Integer
    SobolSampling(n) = n > 0 ? new(n) : error("n must be greater than zero")
end

const Φ = Normal()

function sample(inputs::Array{<:UQInput}, sim::MonteCarlo)
    sample(inputs, sim.n)
end

function sample(inputs::Array{<:UQInput}, sim::SobolSampling)

    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)

    n_rv = count_rvs(random_inputs)

    s = SobolSeq(n_rv)

    n_skip = findlast("1", reverse(bitstring(sim.n - 1)))[1] - 1
    skip(s, 2^n_skip)

    u = hcat([next!(s) for i = 1:sim.n]...)'
    samples = DataFrame(quantile.(Φ, u))

    rename!(samples, names(random_inputs))

    if !isempty(deterministic_inputs)
        samples = hcat(samples, sample(deterministic_inputs, sim.n))
    end

    to_physical_space!(inputs, samples)

    return samples
end