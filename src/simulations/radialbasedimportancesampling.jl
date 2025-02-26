"""
  RadialBasedImportanceSampling(n::Integer, β::Real)

  Used to perform radial-based importance sampling with `n` samples and reliability index `β`.
  If no `β`` or `β=0.0` is passed, a [`FORM`](@ref) analysis will automatically be performed to estimate the reliability index.

  # References

[harbitzEfficientSamplingMethod1986](@cite)
"""
mutable struct RadialBasedImportanceSampling <: AbstractSimulation
    n::Integer
    β::Real
end

function RadialBasedImportanceSampling(n::Integer)
    return RadialBasedImportanceSampling(n, 0.0)
end

function sample(inputs::Vector{<:UQInput}, sim::RadialBasedImportanceSampling)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)

    n_rv = count_rvs(random_inputs)
    rv_names = names(random_inputs)

    samples = DataFrame(
        rv_names .=> eachcol(_samples_outside_beta_sphere(sim.β, n_rv, sim.n))
    )

    if !isempty(deterministic_inputs)
        DataFrames.hcat!(samples, sample(deterministic_inputs, sim.n))
    end

    to_physical_space!(inputs, samples)

    return samples
end

function _samples_outside_beta_sphere(β::Real, k::Integer, n::Integer)

    # sample chi square distributed radius > β
    r = sqrt.(rand(truncated(Chisq(k), β^2, Inf), n))

    # standard normal samples
    ϕ = rand(Normal(), (n, k))
    # normalize each row to create samples on the surface of an n-sphere with radius 1.0
    ϕ = ϕ ./ norm.(eachrow(ϕ))

    # scale samples with the radius to > β
    x = ϕ .* r

    # error if samples are generated inside the sphere
    @assert all(norm.(eachrow(x)) .> β)

    return x
end

Χ
