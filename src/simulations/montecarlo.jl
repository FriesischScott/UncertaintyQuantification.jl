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

struct LatinHypercubeSampling <: AbstractQuasiMonteCarlo
    n::Integer
    LatinHypercubeSampling(n) = n > 0 ? new(n) : error("n must be greater than zero")
end

struct LatticeRuleSampling <: AbstractQuasiMonteCarlo
    n::Integer
    LatticeRuleSampling(n) = n > 0 ? new(n) : error("n must be greater than zero")
end

function sample(inputs::Vector{<:UQInput}, sim::MonteCarlo)
    return sample(inputs, sim.n)
end

function sample(inputs::Vector{<:UQInput}, sim::AbstractQuasiMonteCarlo)
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
    return transpose(QuasiMonteCarlo.sample(sim.n, zeros(rvs), ones(rvs), SobolSample()))
end

function qmc_samples(sim::HaltonSampling, rvs::Integer)
    coprimes = [prime(i) for i in 1:rvs]
    return transpose(
        QuasiMonteCarlo.sample(
            sim.n, zeros(rvs), ones(rvs), LowDiscrepancySample(coprimes, false)
        ),
    )
end

function qmc_samples(sim::LatinHypercubeSampling, rvs::Integer)
    return transpose(
        QuasiMonteCarlo.sample(sim.n, zeros(rvs), ones(rvs), LatinHypercubeSample())
    )
end

function qmc_samples(sim::LatticeRuleSampling, rvs::Integer)
    return transpose(
        QuasiMonteCarlo.sample(sim.n, zeros(rvs), ones(rvs), LatticeRuleSample())
    )
end
