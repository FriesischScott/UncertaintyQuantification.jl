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