function gradient(
    models::Array{<:UQModel},
    inputs::Array{<:UQInput},
    x::DataFrame,
    output::Symbol,
)

    samples = copy(x)

    random_names = _random_inputs(inputs)

    function f(x)
        samples = sample(inputs, size(x, 1))
        samples[:, random_names] .= x

        for m in models
            samples = evaluate(m, samples)
        end

        return samples[:, output][1]
    end

    g = grad(central_fdm(2, 1), f, x)
    names!(g, random_names)

    return g
end

function gradient_in_standard_normal_space(
    models::Array{<:UQModel},
    inputs::Array{<:UQInput},
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

    g = DataFrame(grad(forward_fdm(2, 1), f, reference))
    names!(g, random_names)

    return g
end
