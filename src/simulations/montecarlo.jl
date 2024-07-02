struct MonteCarlo <: AbstractMonteCarlo
    n::Integer
    MonteCarlo(n) = n > 0 ? new(n) : error("n must be greater than zero")
end

struct SobolSampling <: AbstractQuasiMonteCarlo
    n::Integer
    randomization::Symbol

    function SobolSampling(n::Integer, randomization::Symbol=:matousekscramble)
        randomization ∉ [:matousekscramble, :owenscramble, :none] &&
            error("type must be :matousekscramble :owenscramble or :none")
        if n > 0
            if !isinteger(log2(n))
                n = Int(2^ceil(log2(n)))
                @warn("n must be a power of 2, automatically increased to $n")
            end
            return new(n, randomization)
        else
            error("n must be greater than zero")
        end
    end
end

struct FaureSampling <: AbstractQuasiMonteCarlo
    n::Integer
    randomization::Symbol

    function FaureSampling(n::Integer, randomization::Symbol=:matousekScramble)
        randomization ∉ [:matousekscramble, :owenscramble, :none] &&
            error("type must be :matousekscramble, :owenscramble or :none")
        if n > 0
            return new(n, randomization)
        else
            error("n must be greater than zero")
        end
    end
end

struct HaltonSampling <: AbstractQuasiMonteCarlo
    n::Integer
    randomization::Symbol

    function HaltonSampling(n::Integer, randomization::Symbol=:none)
        randomization ∉ [:none] && error("type must be :none")
        if n > 0
            return new(n, randomization)
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

    function LatticeRuleSampling(n::Integer, randomization::Symbol=:shift)
        randomization ∉ [:shift, :none] && error("type must be :shift or :none")
        if n > 0
            return new(n, randomization)
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
        samples = hcat(samples, sample(deterministic_inputs, size(samples, 1)))
    end

    to_physical_space!(inputs, samples)

    return samples
end

sample(input::UQInput, sim::AbstractMonteCarlo) = sample([input], sim)

function qmc_samples(sim::SobolSampling, rvs::Integer)
    return randomize(sim, QuasiMonteCarlo.sample(sim.n, rvs, SobolSample()))
end

function qmc_samples(sim::FaureSampling, rvs::Integer)
    b = nextprime(rvs)
    n = sim.n
    if !isinteger(log(b, sim.n))
        n = Int(b^ceil(log(b, sim.n)))
        @warn(
            "n must be a power of the base (here $b), automatically increased to $n for these samples."
        )
    end
    return randomize(sim, QuasiMonteCarlo.sample(n, rvs, FaureSample()), b)
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

function randomize(sim::AbstractQuasiMonteCarlo, u::Matrix, b=2)
    if sim.randomization == :matousekscramble
        u = QuasiMonteCarlo.randomize(u, MatousekScramble(; base=b))
    elseif sim.randomization == :owenscramble
        u = QuasiMonteCarlo.randomize(u, OwenScramble(; base=b))
    elseif sim.randomization == :shift
        u = QuasiMonteCarlo.randomize(u, Shift())
    end

    return u
end

double_samples(sim::MonteCarlo) = MonteCarlo(2 * sim.n)
double_samples(sim::SobolSampling) = SobolSampling(2 * sim.n, sim.randomization)
double_samples(sim::FaureSampling) = FaureSampling(2 * sim.n, sim.randomization)
double_samples(sim::HaltonSampling) = HaltonSampling(2 * sim.n, sim.randomization)
double_samples(sim::LatinHypercubeSampling) = LatinHypercubeSampling(2 * sim.n)
double_samples(sim::LatticeRuleSampling) = LatticeRuleSampling(2 * sim.n, sim.randomization)
