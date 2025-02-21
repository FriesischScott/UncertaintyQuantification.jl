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

#=
L. E. Blumenson, ‘A Derivation of n-Dimensional Spherical Coordinates’, The American Mathematical Monthly, vol. 67, no. 1, pp. 63–66, 1960, doi: 10.2307/2308932.
=#
function _samples_outside_beta_sphere(β::Real, k::Integer, n::Integer)
    r = sqrt.(rand(truncated(Chisq(k), β^2, Inf), n))
    if k == 1
        return r
    elseif k == 2
        φ = rand(Uniform(0, 2π), n)

        return [r .* cos.(φ) r .* sin.(φ)]
    else
        φ = rand(Uniform(0, π), (n, k - 2))

        θ = rand(Uniform(0, 2π), n)

        x = r .* cos.(φ[:, 1])

        for j in 2:(k - 2)
            x = hcat(x, r .* cos.(φ[:, j]) .* prod(sin.(φ[:, 1:(j - 1)]); dims=2))
        end

        x = hcat(x, r .* sin.(θ) .* prod(sin.(φ[:, 1:(k - 2)]); dims=2))

        x = hcat(x, r .* cos.(θ) .* prod(sin.(φ[:, 1:(k - 2)]); dims=2))

        return x
    end
end
