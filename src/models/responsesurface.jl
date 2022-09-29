using DynamicPolynomials

struct ResponseSurface <: UQModel
    deg::Int64
    dim::Int64
    parameters::Array
    y::Symbol

    function ResponseSurface(data::DataFrame, dependendVarName::Symbol, deg::Int64, dim::Int64)
        if deg < 0 || dim <= 0
            print("degree must be non negative, dimension must be positive")
            return nothing
        end
        
        paramCount = findNumberOfParameters(dim, deg)
        

        if size(data, 1) < paramCount
            print("sample size too small.")
            return nothing
        end

        params = multiDimensionalPolynomialRegression(data, dependendVarName, deg, dim)

        return new(deg, dim, params, dependendVarName)
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

    println(inv(M' * M) * M' * Y)
    return inv(M' * M) * M' * Y
end

function evaluate!(rs::ResponseSurface, data::DataFrame)
    if rs.dim - 1 != size(data,2)
        println("Data does not fit given ResponseSurface")
        return nothing
    end

    data[!, rs.y] = zeros(size(data, 1))
    @polyvar vars[1:rs.dim - 1]
    Mon = reverse(monomials(vars, 0:rs.deg))
    row = ones(1, size(Mon, 1))

    for i in 1:size(data, 1)
        sub_mon = Mon
        if i % 50 == 0
            print("1111------------")
            println(sub_mon)
        end
        for j in 1:rs.dim - 1
            sub_mon = subs(sub_mon, vars[j]=>data[i, j])
        end
        if i % 50 == 0
            print("2222------------")
            println(sub_mon)
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