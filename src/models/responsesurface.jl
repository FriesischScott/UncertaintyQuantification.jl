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

julia> rs = ResponseSurface(data, :y, 2, 2)
ResponseSurface([0.4833333333331211, -0.23863636363637397, 1.018939393939391], :y, 2, 2)
```
"""
struct ResponseSurface <: UQModel
    parameters::Array
    y::Symbol
    deg::Int64
    dim::Int64

    function ResponseSurface(data::DataFrame, dependendVarName::Symbol, deg::Int64, dim::Int64)
        if deg < 0 || dim <= 0
            print("degree must be non negative, dimension must be positive")
            return nothing
        end

        if size(data, 1) < findNumberOfParameters(dim, deg)
            print("sample size too small.")
            return nothing
        end

        params = multiDimensionalPolynomialRegression(data, dependendVarName, deg, dim)

        return new(params, dependendVarName, deg, dim)
    end
end



#only to be called internally
function findNumberOfParameters(dim::Int64, deg::Int64)
    if deg == 0
        return 1
    else 
        # mit Wiederholung, ohne Sortierung; n = dim -1; k = deg
        return Int((factorial(dim + deg -2) / (factorial(deg) * factorial(dim - 2))) + findNumberOfParameters(dim, deg - 1)) 
    end
end



# only to be called internally by constructor
function multiDimensionalPolynomialRegression(data::DataFrame, y::Symbol, deg::Int64, dim::Int64)
    if size(data, 2) != dim
        println("Dataframe does not fit given Dimension.")
        return nothing
    end

    Y = data[:, y]
    x_data = (data[:, Not(String(y))])

    @polyvar vars[1:dim - 1]
    Mon = reverse(monomials(vars, 0:deg))
    row = Mon
    # matrix for final computation
    M = ones(size(data, 1), (size(Mon, 1)))


    # manually doing 1st iteration of following loop, so that X is defined
    for j in 1:(dim - 1)
        row = subs(row, vars[j]=>x_data[1, j])
    end
    X = row'

    #filling X with polynomials and inserting data
    for i in 2:size(data, 1)
        row = Mon

        for j in 1:(dim - 1)
            row = subs(row, vars[j]=>x_data[i, j])
        end

        X = vcat(X, row')

    end

    #copying X into M, to change entry type from polyvar to Float; also adding column of ones
    for i in 1:size(X, 1)
        for j in 1:size(X, 2)
            M[i, j] = X[i, j]
        end
    end

    return inv(M' * M) * M' * Y
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

#output

3×2 DataFrame
 Row │ x        y
     │ Float64  Float64   
─────┼────────────────────
   1 │     2.5    6.25511
   2 │    11.0  121.15
   3 │    15.0  226.165
```
"""
function evaluate!(rs::ResponseSurface, data::DataFrame)
    if rs.dim - 1 != size(data,2)
        println("Data does not fit given ResponseSurface")
        return nothing
    end

    data[!, rs.y] = zeros(size(data, 1))    #result column

    
    @polyvar vars[1:rs.dim - 1]
    Mon = reverse(monomials(vars, 0:rs.deg))
    row = ones(1, size(Mon, 1))

    for i in 1:size(data, 1)
        sub_mon = Mon   #"copying" monomials

        for j in 1:rs.dim - 1
            sub_mon = subs(sub_mon, vars[j]=>data[i, j])    #filling 1 row of monomials with data
        end

        #typeconverting the data-filled monomials; transposing the vector
        for k in 1:size(row, 2)
            row[1, k] = sub_mon[k, 1] 
        end
       
        val = row * rs.parameters
        data[i, rs.y] = val[1, 1]
    end
    return data
end