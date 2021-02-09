struct MonteCarlo <: AbstractMonteCarlo
    n::Integer
    MonteCarlo(n) = n > 0 ? new(n) : error("n must be greater than zero")
end

struct SobolSampling <: AbstractQuasiMonteCarlo
    n::Integer
    SobolSampling(n) = n > 0 ? new(n) : error("n must be greater than zero")
end

struct HaltonSampling <: AbstractQuasiMonteCarlo
    n::Integer
    HaltonSampling(n) = n > 0 ? new(n) : error("n must be greater than zero")
end

function sample(inputs::Array{<:UQInput}, sim::MonteCarlo)
    sample(inputs, sim.n)
end

function sample(inputs::Array{<:UQInput}, sim::AbstractQuasiMonteCarlo)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)

    n_rv = count_rvs(random_inputs)

    u = qmc_samples(sim, n_rv)

    samples = quantile.(Normal(), u)
    samples = DataFrame(names(random_inputs) .=> eachcol(samples))

    if !isempty(deterministic_inputs)
        samples = hcat(samples, sample(deterministic_inputs, sim.n))
    end

    to_physical_space!(inputs, samples)

    return samples
end

function qmc_samples(sim::SobolSampling, rvs::Integer)
    s = SobolSeq(rvs)

    n_skip = findlast("1", reverse(bitstring(sim.n - 1)))[1] - 1
    skip(s, 2^n_skip)

    u = hcat([next!(s) for i = 1:sim.n]...) |> transpose
end

function qmc_samples(sim::HaltonSampling, rvs::Integer)
    h = HaltonPoint(rvs, length=sim.n)

    u = zeros(sim.n, rvs)

    for (i, hp) ∈ enumerate(h)
        u[i, :] = hp
    end

    return u
end