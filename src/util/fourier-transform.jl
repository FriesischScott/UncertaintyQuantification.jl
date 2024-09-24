function dft(x::AbstractVector)
    N = length(x)
    X = zeros(ComplexF64, N)

    for n in 1:N
        X += x[n] * exp.(-1im * 2 * pi .* (0:N-1) ./ N * (n-1))
    end

    return X
end

function idft(X::AbstractVector)
    N = length(X)
    x = zeros(ComplexF64, N)

    for k in 1:N
        x += 1/N * X[k] * exp.(1im * 2 * pi * (k-1) / N .* (0:N-1))
    end

    return x
end
