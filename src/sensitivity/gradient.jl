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
