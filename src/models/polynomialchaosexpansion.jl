struct PolynomialChaosExpansion <: UQModel
    S::Array{Float64}
    n::Array{Symbol}
    P::Integer
    inputs::Array{UQInput,1}
    output::Symbol

    function PolynomialChaosExpansion(
        data::DataFrame,
        inputs::Array{UQInput,1},
        P::Integer,
        output::Symbol
        )

        random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
        random_names = names(random_inputs)

        x = data[:, random_names]
        to_standard_normal_space!(random_inputs, x)

        N = size(x, 1)
        M = size(x, 2)

        A = zeros(size(x, 1), P)

        for i ∈ 1:N, j ∈ 1:P
            A[i, j] = Ψ(convert(Array{Float64,1}, x[i, random_names]), j - 1)
        end

        S = inv(transpose(A) * A) * transpose(A) * data[:, output]

        new(S, random_names, P, random_inputs, output)
    end
end

function evaluate!(
    pce::PolynomialChaosExpansion,
    df::DataFrame
    )

    data = df[:, pce.n]
    to_standard_normal_space!(pce.inputs, data)
    # convert to matrix and order variables by ps.n
    x = convert(Matrix, data)

    out = map(row -> dot(pce.S, [Ψ(convert(Matrix, row), j) for j ∈ 0:pce.P - 1]), eachrow(x))
    df[!, pce.output] = out
end

function Ψ(x::Array{Float64,1}, n::Integer)
    d = repeat([n], length(x))
    ops = GaussOrthoPoly.(d)
    mop = MultiOrthoPoly(ops, n)

    PolyChaos.evaluate(x, mop) |> sum
end