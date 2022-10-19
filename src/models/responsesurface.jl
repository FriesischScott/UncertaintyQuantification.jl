using DynamicPolynomials
"""
    ResponseSurface(data::DataFrame, dependendVarName::Symbol, deg::Int64, dim::Int64)

Creates a response surface using polynomial least squares regression with given degree and dimension

# Examples
```jldoctest
julia> data = DataFrame(x = 1:10, y = [1, 4, 10, 15, 24, 37, 50, 62, 80, 101])
10×2 DataFrame
 Row │ x      y     
     │ Int64  Int64 
─────┼──────────────
   1 │     1      1
   2 │     2      4
   3 │     3     10
   4 │     4     15
   5 │     5     24
   6 │     6     37
   7 │     7     50
   8 │     8     62
   9 │     9     80
  10 │    10    101

julia> rs = ResponseSurface(data, :y, 2)
ResponseSurface([1.018939393939398, -0.23863636363631713, 0.4833333333332348], :y, [:x], 2, Monomial{true}[x₁², x₁, 1])
```
"""
struct ResponseSurface <: UQModel
    parameters::Array
    y::Symbol
    names::Array{Symbol}
    deg::Int64
    monomial::DynamicPolynomials.MonomialVector{true}

    function ResponseSurface(data::DataFrame, dependendVarName::Symbol, deg::Int64)
        if deg < 0
            print("degree must be non negative")
            return nothing
        end

        @polyvar x[1:size(data, 2) - 1]
        monomial = monomials(x, 0:deg)

        names = propertynames(data[:, Not(dependendVarName)])

        params = multiDimensionalPolynomialRegression(data, dependendVarName, monomial)

        return new(params, dependendVarName, names, deg, monomial)
    end
end



# only to be called internally by constructor
function multiDimensionalPolynomialRegression(data::DataFrame, y::Symbol, monomial::DynamicPolynomials.MonomialVector{true})
    Y = data[:, y]
    x_data = Matrix{Float64}(data[:, Not(String(y))]) # convert to matrix
    M = zeros(size(data, 1), size(monomial, 1))#prepare matrix

    for i in 1 : size(M, 1)
        M[i, :] = map(m -> m(convert(Array, x_data[i, :])), monomial') #fill matrix with monomials filled with given data
    end

    return inv(transpose(M) * M) * transpose(M) * Y #regression formula
end



#called internally by evaluate!
#evaluates one datapoint using a given ResponseSurface
function calc(row::Array, rs::ResponseSurface)
    map(m -> m(row), rs.monomial') * rs.parameters
end



"""
evaluate!(rs::ResponseSurface, data::DataFrame)
    evaluating data by using a previously trained ResponseSurface


#Examples

```jldoctest
julia> data = DataFrame(x = 1:10, y = [1, 4, 10, 15, 24, 37, 50, 62, 80, 101])
10×2 DataFrame
 Row │ x      y     
     │ Int64  Int64 
─────┼──────────────
   1 │     1      1
   2 │     2      4
   3 │     3     10
   4 │     4     15
   5 │     5     24
   6 │     6     37
   7 │     7     50
   8 │     8     62
   9 │     9     80
  10 │    10    101

julia> rs = ResponseSurface(data, :y, 2, 2)
ResponseSurface([0.4833333333331211, -0.23863636363637397, 1.018939393939391], :y, 2, 2)
julia> evaluate!(rs, [2.5, 11, 15])
```
"""
function evaluate!(rs::ResponseSurface, data::DataFrame)
    x = Matrix{Float64}(data[:, rs.names]) # convert to matrix, sort by rs.names
    out = map(row -> (calc(convert(Array, row), rs)), eachrow(x)) # fill monomial, evaluate given data with ResponseSurface
    data[!, rs.y] = out
    return nothing
end

