function gradient(
    models::Array{<:AbstractModel},
    inputs::Array{<:AbstractInput},
    x,
    output::Symbol,
)

    random_names = Vector{Symbol}()

    for i in inputs
        if isa(i, RandomVariable)
            push!(random_names, Symbol(i.name))
        elseif isa(i, RandomVariableSet)
            append!(random_names, map(x -> Symbol(x.name), i.members))
        end
    end

    function f(x)
        samples = sample(inputs, size(x, 1))
        samples[:, random_names] .= x

        for m in models
            samples = evaluate(m, samples)
        end

        return samples[:, output][1]
    end

    grad(central_fdm(2, 1), f, x)
end

function gradient_in_standard_normal_space(
    models::Array{<:AbstractModel},
    inputs::Array{<:AbstractInput},
    x::DataFrame,
    output::Symbol,
)

    samples = copy(x)

    random_names = _random_inputs(inputs)
    to_standard_normal_space!(inputs, samples)

    function f(x)
        samples[:, random_names] .= x
        to_physical_space!(inputs, samples)
        for m in models
            samples = evaluate(m, samples)
        end
        return samples[:, output][1]
    end

    reference = convert(Matrix, samples[:, random_names])

    g = DataFrame(-1 * grad(forward_fdm(2, 1), f, reference))
    names!(g, random_names)

    return g
end
