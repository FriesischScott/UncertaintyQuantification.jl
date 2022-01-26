struct PolynomialChaosExpansion <: UQModel
    y::Vector{Float64}
    n::Vector{Symbol}
    p::Integer
    inputs::Vector{<:UQInput}
    output::Symbol

    function PolynomialChaosExpansion(
        data::DataFrame, inputs::Vector{<:UQInput}, p::Integer, output::Symbol
    )
        random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
        random_names = names(random_inputs)

        x = data[:, random_names]
        to_standard_normal_space!(random_inputs, x)

        N = size(x, 1)

        P = Int(factorial(p + size(x, 2)) / (factorial(p) * factorial(size(x, 2))))

        A = zeros(size(x, 1), P)

        for i in 1:N, j in 1:P
            A[i, j] = Ψ(collect(x[i, :]), j - 1)
        end

        y = inv(transpose(A) * A) * transpose(A) * data[:, output]

        return new(y, random_names, P, random_inputs, output)
    end
end

function evaluate!(pce::PolynomialChaosExpansion, df::DataFrame)
    data = df[:, pce.n]
    to_standard_normal_space!(pce.inputs, data)

    p = [0:(length(pce.y) - 1);]
    out = map(row -> sum(pce.y .* [Ψ(collect(row), j) for j in p]), eachrow(data))
    return df[!, pce.output] = out
end

function Ψ(x::Vector{Float64}, n::Integer)
    if n == 0
        return 1
    end

    d = repeat([n], length(x))
    ops = GaussOrthoPoly.(d)
    mop = MultiOrthoPoly(ops, n)

    return sum(PolyChaos.evaluate(x, mop))
end
