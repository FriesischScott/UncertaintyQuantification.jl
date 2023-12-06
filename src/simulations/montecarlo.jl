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

function sample(inputs::Vector{<:PreciseUQInput}, sim::MonteCarlo)
    return sample(inputs, sim.n)
end

function sample(inputs::Vector{<:PreciseUQInput}, sim::AbstractQuasiMonteCarlo)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)

    n_rv = count_rvs(random_inputs)

    u = qmc_samples(sim, n_rv)

    samples = quantile.(Normal(), u)
    samples = DataFrame(names(random_inputs) .=> eachrow(samples))

    if !isempty(deterministic_inputs)
        samples = hcat(samples, sample(deterministic_inputs, sim.n))
    end

    to_physical_space!(inputs, samples)

    return samples
end

sample(input::PreciseUQInput, sim::AbstractMonteCarlo) = sample([input], sim)

function qmc_samples(sim::SobolSampling, rvs::Integer)
    return QuasiMonteCarlo.sample(sim.n, rvs, SobolSample())
end

function qmc_samples(sim::HaltonSampling, rvs::Integer)
    samples = QuasiMonteCarlo.sample(sim.n, rvs, HaltonSample())
    return rvs > 1 ? samples : reshape(samples, 1, sim.n)
end

function qmc_samples(sim::LatinHypercubeSampling, rvs::Integer)
    return QuasiMonteCarlo.sample(sim.n, rvs, LatinHypercubeSample())
end

function qmc_samples(sim::LatticeRuleSampling, rvs::Integer)
    return QuasiMonteCarlo.sample(sim.n, rvs, LatticeRuleSample())
end

double_samples(sim::MonteCarlo) = MonteCarlo(2 * sim.n)
double_samples(sim::SobolSampling) = SobolSampling(2 * sim.n)
double_samples(sim::HaltonSampling) = HaltonSampling(2 * sim.n)
double_samples(sim::LatinHypercubeSampling) = LatinHypercubeSampling(2 * sim.n)
double_samples(sim::LatticeRuleSampling) = LatticeRuleSampling(2 * sim.n)
