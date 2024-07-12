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
    performance::Function;
    fdm::FiniteDifferencesMethod=CentralFiniteDifferences(3),
)
    samples = copy(reference)

    random_names = names(
        filter(i -> isa(i, RandomUQInput) || isa(i, ProbabilityBox), inputs)
    )

    function f(x)
        samples_eval = deepcopy(reference)
        samples_eval[:, random_names] .= x

        to_physical_space!(inputs, samples_eval)

        evaluate!(models, samples_eval)

        return performance(samples_eval)[1]
    end

    reference_point = Matrix{Float64}(samples[:, random_names])

    g = _grad(f, fdm, reference_point)
    return (; zip(random_names, g)...)
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
