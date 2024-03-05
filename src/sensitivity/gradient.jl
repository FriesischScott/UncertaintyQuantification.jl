function gradient(
    models::Vector{<:UQModel},
    inputs::Vector{<:UQInput},
    x::DataFrame,
    output::Symbol;
    fdm::FiniteDifferencesMethod=CentralFiniteDifferences(3),
)
    samples = copy(x)

    random_names = names(filter(i -> isa(i, RandomUQInput), inputs))

    function f(x)
        samples = sample(inputs, size(x, 1))
        samples[:, random_names] .= x

        evaluate!(models, samples)

        return samples[:, output][1]
    end

    reference = Matrix{Float64}(samples[:, random_names])

    g = _grad(f, fdm, reference)

    return (; zip(random_names, g)...)
end

function gradient_in_standard_normal_space(
    models::Vector{<:UQModel},
    inputs::Vector{<:UQInput},
    reference::DataFrame,
    output::Symbol;
    fdm::FiniteDifferencesMethod=CentralFiniteDifferences(3),
)
    samples = copy(reference)

    random_names = names(filter(i -> isa(i, RandomUQInput), inputs))
    to_standard_normal_space!(inputs, samples)

    function f(x)
        samples[:, random_names] .= x
        to_physical_space!(inputs, samples)

        evaluate!(models, samples)

        return samples[:, output][1]
    end

    reference = Matrix{Float64}(samples[:, random_names])

    g = _grad(f, fdm, reference)

    return (; zip(random_names, g)...)
end

function gradient_in_metric_space(
    models::Vector{<:UQModel},
    inputs::Vector{<:PreciseUQInput},
    reference::DataFrame,
    output::Symbol,
    metric::Function;
    fdm::FiniteDifferencesMethod=CentralFiniteDifferences(3),
)
    samples = copy(reference)

    random_names = names(filter(i -> isa(i, RandomUQInput), inputs))
    to_standard_normal_space!(inputs, samples)



    function f(x)
        x_transformed = metric_transformation(x, metric, LinearAlgebra.norm)
        samples[:, random_names] .= x_transformed
        to_physical_space!(inputs, samples)

        evaluate!(models, samples)

        return samples[:, output][1]
    end

    reference = Matrix{Float64}(samples[:, random_names])

    reference_metric = metric_transformation(reference, LinearAlgebra.norm, metric)

    g = _grad(f, fdm, reference_metric)

    return (; zip(random_names, g)...)    

end

function metric_transformation(samples, norm1, norm2)
    if iszero(samples) return samples end

    Ns = size(samples, 1)
    norm_p1 = [norm1(samples[i,:])  for i in 1:Ns] 
    norm_p1 = reduce(vcat, norm_p1)

    samples_normlised = samples ./ norm_p1

    norm_p2 = [norm2(samples[i,:])  for i in 1:Ns] 

    return samples_normlised .* norm_p2
end

function _grad(f::Function, fdm::CentralFiniteDifferences, x)
    return grad(central_fdm(fdm.order, fdm.derivative), f, x)[1]
end

function _grad(f::Function, fdm::ForwardFiniteDifferences, x)
    return grad(forward_fdm(fdm.order, fdm.derivative), f, x)[1]
end

function _grad(f::Function, fdm::BackwardFiniteDifferences, x)
    return grad(backward_fdm(fdm.order, fdm.derivative), f, x)[1]
end
