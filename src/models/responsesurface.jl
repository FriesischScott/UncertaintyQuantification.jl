"""
    ResponseSurface(data::DataFrame, dependendVarName::Symbol, deg::Int, dim::Int)

Creates a response surface using polynomial least squares regression with given degree.

# Examples
```jldoctest
julia> data = DataFrame(x = 1:10, y = [1, 4, 10, 15, 24, 37, 50, 62, 80, 101]);

julia> rs = ResponseSurface(data, :y, 2)
ResponseSurface([0.48333333333332457, -0.23863636363636026, 1.0189393939393936], :y, [:x], 2, Monomials.Monomial[1, x1, x1²])
```
"""
struct ResponseSurface <: UQModel
    β::Array
    y::Symbol
    names::Vector{Symbol}
    p::Int64
    monomials::Vector{Monomial}

    function ResponseSurface(data::DataFrame, output::Symbol, p::Int)
        if p < 0
            error("Degree(p) of ResponseSurface must be non-negative.")
        end

        x = ["x$i" for i in 1:(size(data, 2) - 1)]
        m = monomials(x, p, GradedLexicographicOrder(); include_zero=true)

        names = propertynames(data[:, Not(output)])

        X = transpose(Matrix{Float64}(data[:, names])) # convert to matrix, sort by rs.names
        y = Vector{Float64}(data[:, output])

        β = m(X)' \ y
        return new(β, output, names, p, m)
    end
end

"""
    evaluate!(rs::ResponseSurface, data::DataFrame)

evaluating data by using a previously trained ResponseSurface.


# Examples

```jldoctest
julia> data = DataFrame(x = 1:10, y = [1, 4, 10, 15, 24, 37, 50, 62, 80, 101]);

julia> rs = ResponseSurface(data, :y, 2);

julia> df = DataFrame(x = [2.5, 11, 15]);

julia> evaluate!(rs, df);

julia> df.y |> DisplayAs.withcontext(:compact => true)
3-element Vector{Float64}:
   6.25511
 121.15
 226.165
```
"""
function evaluate!(rs::ResponseSurface, data::DataFrame)
    x = transpose(Matrix{Float64}(data[:, rs.names])) # convert to matrix, sort by rs.names
    out = map(x -> dot(x, rs.β), eachcol(rs.monomials(x)))

    data[!, rs.y] = out
    return nothing
end
