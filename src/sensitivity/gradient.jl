function gradient(
    models::Array{<:UQModel}, inputs::Array{<:UQInput}, x::DataFrame, output::Symbol
)
    samples = copy(x)

    random_names = names(filter(i -> isa(i, RandomUQInput), inputs))

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
    models::Array{<:UQModel}, inputs::Array{<:UQInput}, x::DataFrame, output::Symbol
)
    random_names = names(filter(i -> isa(i, RandomUQInput), inputs))

    function f(x)
        sample = DataFrame(random_names .=> x)
        to_physical_space!(inputs, sample)

        evaluate!(models, sample)

        return sample[1, output][1]
    end

    reference = copy(x)
    to_standard_normal_space!(inputs, reference)
    reference = vec(Matrix(reference[:, random_names]))

    g = grad(forward_fdm(2, 1), f, reference)[1]

    return (; zip(random_names, g)...)
end
