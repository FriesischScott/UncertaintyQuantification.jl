"""
    ResponseSurface(data::DataFrame, dependendVarName::Symbol, deg::Int, dim::Int)

Creates a response surface using polynomial least squares regression with given degree.

# Examples
```jldoctest
julia> data = DataFrame(x = 1:10, y = [1, 4, 10, 15, 24, 37, 50, 62, 80, 101]);

julia> rs = ResponseSurface(data, :y, 2) |> DisplayAs.withcontext(:compact => true)
ResponseSurface([1.01894, -0.238636, 0.483333], :y, [:x], 2, DynamicPolynomials.Monomial{true}[x₁², x₁, 1])
```
"""
struct ResponseSurface <: UQModel
    β::Array
    y::Symbol
    names::Vector{Symbol}
    p::Int
    monomials::MonomialVector{true}

    function ResponseSurface(data::DataFrame, output::Symbol, p::Int)
        if p < 0
            error("Degree(p) of ResponseSurface must be non-negative.")
        end

        @polyvar x[1:(size(data, 2) - 1)]
        m = monomials(x, 0:p)

        names = propertynames(data[:, Not(output)])

        X = Matrix{Float64}(data[:, names]) # convert to matrix, sort by rs.names
        y = Vector{Float64}(data[:, output])

        β = multi_dimensional_polynomial_regression(X, y, m)

        return new(β, output, names, p, m)
    end
end

# only to be called internally by constructor
function multi_dimensional_polynomial_regression(
    X::Matrix, y::Vector, monomials::MonomialVector{true}
)

    #fill monomials with the given x values for each row
    M = mapreduce(row -> begin
        return map(m -> m(row), monomials)'
    end, vcat, eachrow(X))

    return M \ y   #β
end

#called internally by evaluate!
#evaluates one datapoint using a given ResponseSurface
function calc(row::Array, rs::ResponseSurface)
    return map(m -> m(row), rs.monomials') * rs.β
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
    x = Matrix{Float64}(data[:, rs.names]) # convert to matrix, sort by rs.names
    out = map(row -> (calc(convert(Array, row), rs)), eachrow(x)) # fill monomial, evaluate given data with ResponseSurface
    data[!, rs.y] = out
    return nothing
end
