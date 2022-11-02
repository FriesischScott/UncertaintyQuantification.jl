"""
    PolyharmonicSpline(data::DataFrame, k::Int64, output::Symbol)

creates a polyharmonic spline that is trained by given data

#Examples
```jldoctest

julia> data = DataFrame(x = 1:10, y = [1, -5, -10, -12, -8, -1, 5, 12, 23, 50])
10×2 DataFrame
 Row │ x      y
     │ Int64  Int64
─────┼──────────────
   1 │     1      1
   2 │     2     -5
   3 │     3    -10
   4 │     4    -12
   5 │     5     -8
   6 │     6     -1
   7 │     7      5
   8 │     8     12
   9 │     9     23
  10 │    10     50

  julia> PolyharmonicSpline(data, 2, :y)
PolyharmonicSpline([1.1473268119780278; -0.44960947031466086; … ; -5.331010776968267; 3.8862763174093313;;],
     [-112.00528786482354; 6.844431546357826;;], [1.0; 2.0; … ; 9.0; 10.0;;], 2, [:x], :y)
```
"""
struct PolyharmonicSpline <: UQModel
    w::Array{Float64}
    v::Array{Float64}
    c::Array{Float64,2}
    k::Int64
    n::Array{Symbol}
    output::Symbol

    function PolyharmonicSpline(data::DataFrame, k::Int64, output::Symbol)
        f = data[:, output]

        centers = select(data, Not(output))
        n = propertynames(centers)
        centers = Matrix{Float64}(centers)

        dim = size(centers, 1)

        A = zeros(dim, dim)
        for i in 1:dim, j in 1:dim
            if i == j
                continue
            end
            r = sqrt(sum((centers[i, :] - centers[j, :]) .^ 2))
            A[i, j] = ϕ(r, k)
        end

        B = [ones(dim, 1) centers]

        M = [A B; transpose(B) zeros(size(B, 2), size(B, 2))]
        F = [f; zeros(size(B, 2), 1)]

        wv = M \ F

        w = wv[1:dim, :]
        v = wv[(dim + 1):size(wv, 1), :]

        return new(w, v, centers, k, n, output)
    end
end

function ϕ(r::Float64, k::Int64)
    if k % 2 != 0
        return r^k
    elseif r < 1
        return r^(k - 1) * log(r^r)
    else
        return r^k * log(r)
    end
end

function calc(ps::PolyharmonicSpline, x::Array{Float64,1})
    r = sqrt.(sum((ps.c .- transpose(x)) .^ 2; dims=2))
    f = sum(ϕ.(r, ps.k) .* ps.w)
    return f += (transpose(ps.v) * [1; x])[1]
end

"""
    evaluate!(ps::PolyharmonicSpline, df::DataFrame)

    evaluates given data using a previously contructed PolyharmonicSpline

#Examples
```jldoctest

data = DataFrame(x = 1:10, y = [1, -5, -10, -12, -8, -1, 5, 12, 23, 50])

ps = PolyharmonicSpline(data, 2, :y)

df = DataFrame( x = [2.5, 7.5, 12, 30])

evaluate!(ps, df)

# output

4-element Vector{Float64}:
  -7.754272281066534
   8.290831024829075
  84.4685159898265
 260.4367316915062
```
"""
function evaluate!(ps::PolyharmonicSpline, df::DataFrame)
    x = Matrix{Float64}(df[:, ps.n]) # convert to matrix and order variables by ps.n

    out = map(row -> calc(ps, convert(Array, row)), eachrow(x))
    return df[!, ps.output] = out
end
