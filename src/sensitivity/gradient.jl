function gradient(
    models::Array{<:UQModel},
    inputs::Array{<:UQInput},
    x::DataFrame,
    output::Symbol
)
    samples = copy(x)

    random_names = filter(i -> isa(i, RandomUQInput), inputs) |> names

    function f(x)
        samples = sample(inputs, size(x, 1))
        samples[:, random_names] .= x

        for m in models
            evaluate!(m, samples)
        end

        return samples[:, output][1]
    end

    reference = convert(Array{Float64, 2}, samples[:, random_names])

    g = DataFrame(grad(forward_fdm(2, 1), f, reference)[1])
    rename!(g, random_names)

    return g
end

function gradient_in_standard_normal_space(
    models::Array{<:UQModel},
    inputs::Array{<:UQInput},
    reference::DataFrame,
    output::Symbol
)
    samples = copy(reference)

    random_names = filter(i -> isa(i, RandomUQInput), inputs) |> names
    to_standard_normal_space!(inputs, samples)

    function f(x)
        samples[:, random_names] .= x
        to_physical_space!(inputs, samples)
        for m in models
            evaluate!(m, samples)
        end
        return samples[:, output][1]
    end

    reference = convert(Array{Float64, 2}, samples[:, random_names])

    g = DataFrame(grad(forward_fdm(2, 1), f, reference)[1])
    rename!(g, random_names)

    return g
end
