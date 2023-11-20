struct MonteCarlo <: AbstractMonteCarlo
    n::Integer
    MonteCarlo(n) = n > 0 ? new(n) : error("n must be greater than zero")
end

struct SobolSampling <: AbstractQuasiMonteCarlo
    n::Integer
    randomization::Symbol
    base::Integer
    pad::Integer

    function SobolSampling(n::Integer, randomization::Symbol=:matousekscramble)
        randomization ∉ [:matousekscramble, :digitalshift, :shift, :owenscramble, :none] &&
            error(
                "type must be :matousekscramble, :digitalshift, :shift, :owenscramble or :none",
            )
        return if n > 0
            new(n, randomization, 2, 32)
        else
            error("n must be greater than zero")
        end
    end
end

struct HaltonSampling <: AbstractQuasiMonteCarlo
    n::Integer
    randomization::Symbol
    base::Integer
    pad::Integer

    function HaltonSampling(
        n::Integer,
        randomization::Symbol=:matousekscramble,
        base::Integer=2,
        pad::Integer=32,
    )
        randomization ∉ [:matousekscramble, :digitalshift, :shift, :owenscramble, :none] &&
            error(
                "type must be :matousekscramble, :digitalshift, :shift, :owenscramble or :none",
            )
        pad < log(base, n) && error("pad must be ≥ log(base, n)")
        return if n > 0
            new(n, randomization, base, pad)
        else
            error("n must be greater than zero")
        end
    end
end

struct LatinHypercubeSampling <: AbstractQuasiMonteCarlo
    n::Integer
    LatinHypercubeSampling(n) = n > 0 ? new(n) : error("n must be greater than zero")
end

struct LatticeRuleSampling <: AbstractQuasiMonteCarlo
    n::Integer
    randomization::Symbol
    base::Integer
    pad::Integer

    function LatticeRuleSampling(
        n::Integer,
        randomization::Symbol=:matousekscramble,
        base::Integer=2,
        pad::Integer=32,
    )
        randomization ∉ [:matousekscramble, :digitalshift, :shift, :owenscramble, :none] &&
            error(
                "type must be :matousekscramble, :digitalshift, :shift, :owenscramble or :none",
            )
        pad < log(base, n) && error("pad must be ≥ log(base, n)")
        return if n > 0
            new(n, randomization, base, pad)
        else
            error("n must be greater than zero")
        end
    end
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
    samples = DataFrame(names(random_inputs) .=> eachrow(samples))

    if !isempty(deterministic_inputs)
        samples = hcat(samples, sample(deterministic_inputs, sim.n))
    end

    to_physical_space!(inputs, samples)

    return samples
end

sample(input::UQInput, sim::AbstractMonteCarlo) = sample([input], sim)

function qmc_samples(sim::SobolSampling, rvs::Integer)
    return randomize(sim, QuasiMonteCarlo.sample(sim.n, rvs, SobolSample()))
end

function qmc_samples(sim::HaltonSampling, rvs::Integer)
    samples = QuasiMonteCarlo.sample(sim.n, rvs, HaltonSample())
    return randomize(sim, rvs > 1 ? samples : reshape(samples, 1, sim.n))
end

function qmc_samples(sim::LatinHypercubeSampling, rvs::Integer)
    return QuasiMonteCarlo.sample(sim.n, rvs, LatinHypercubeSample())
end

function qmc_samples(sim::LatticeRuleSampling, rvs::Integer)
    return randomize(sim, QuasiMonteCarlo.sample(sim.n, rvs, LatticeRuleSample()))
end

function randomize(sim::AbstractQuasiMonteCarlo, u::Matrix)
    if sim.randomization == :matousekscramble
        u = QuasiMonteCarlo.randomize(u, MatousekScramble(; base=sim.base, pad=sim.pad))
    elseif sim.randomization == :owenscramble
        u = QuasiMonteCarlo.randomize(u, OwenScramble(; base=sim.base, pad=sim.pad))
    elseif sim.randomization == :digitalshift
        u = QuasiMonteCarlo.randomize(u, DigitalShift(; base=sim.base, pad=sim.pad))
    elseif sim.randomization == :shift
        u = QuasiMonteCarlo.randomize(u, Shift())
    end

    return u
end

double_samples(sim::MonteCarlo) = MonteCarlo(2 * sim.n)
double_samples(sim::SobolSampling) = SobolSampling(2 * sim.n)
double_samples(sim::HaltonSampling) = HaltonSampling(2 * sim.n)
double_samples(sim::LatinHypercubeSampling) = LatinHypercubeSampling(2 * sim.n)
double_samples(sim::LatticeRuleSampling) = LatticeRuleSampling(2 * sim.n)
