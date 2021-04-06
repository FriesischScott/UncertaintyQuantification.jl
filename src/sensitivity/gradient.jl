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

        evaluate!(models, samples)

        return samples[:, output][1]
    end

    reference = Array{Float64,2}(samples[:, random_names])

    g = grad(forward_fdm(2, 1), f, reference)[1]

    return (; zip(random_names, g)...)
end

function gradient_in_standard_normal_space(
    models::Array{<:UQModel},
    inputs::Array{<:UQInput},
    x::DataFrame,
    output::Symbol
)
    samples = copy(x)

    random_names = filter(i -> isa(i, RandomUQInput), inputs) |> names
    to_standard_normal_space!(inputs, samples)

    function f(x)
        samples[:, random_names] .= x
        to_physical_space!(inputs, samples)

        evaluate!(models, samples)

        return samples[:, output][1]
    end

    reference = Array{Float64,2}(samples[:, random_names])

    g = grad(forward_fdm(2, 1), f, reference)[1]

    return (; zip(random_names, g)...)
end
