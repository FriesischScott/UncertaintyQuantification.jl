struct PolynomialChaosExpansion <: UQModel
    y::Vector{Float64}
    n::Vector{Symbol}
    p::Integer
    Ψ::MultiOrthoPoly
    inputs::Vector{<:UQInput}
    output::Symbol

    function PolynomialChaosExpansion(
        data::DataFrame, inputs::Vector{<:UQInput}, p::Integer, output::Symbol
    )
        random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
        random_names = names(random_inputs)

        x = data[:, random_names]
        to_standard_normal_space!(random_inputs, x)

        univariate_polynomials = GaussOrthoPoly.(repeat([p], size(x, 2)))
        Ψ = MultiOrthoPoly(univariate_polynomials, p)

        A = mapreduce(row -> PolyChaos.evaluate(collect(row), Ψ), hcat, eachrow(x))'

        y = inv(transpose(A) * A) * transpose(A) * data[:, output]

        return new(y, random_names, p, Ψ, random_inputs, output)
    end
end

function evaluate!(pce::PolynomialChaosExpansion, df::DataFrame)
    data = df[:, pce.n]
    to_standard_normal_space!(pce.inputs, data)

    out = map(row -> dot(pce.y, PolyChaos.evaluate(collect(row), pce.Ψ)), eachrow(data))
    return df[!, pce.output] = out
end
