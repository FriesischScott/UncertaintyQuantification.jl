struct PolySpline <: UQModel
    w::Array{Float64}
    v::Array{Float64}
    ci::Array{Float64}
    k::Int64
    orderfunc::Function
    name::Symbol
end

function phi(
    x::Vector{},
    c::Vector{},
    k::Int64
    )
    r = sqrt(dot((x .- c), (x .- c)))
    if k%2 != 0
        return r^k
    elseif r < 1
        return r^(k-1) * log(r^r)
    else
        return r^k *log(r)
    end
end

function createpolyspline(
    inputs::Array{},
    f::Array{},
    k::Int64,
    name::Symbol,
    orderfunc::Function = (df -> []) # needed to pull out the correct variables from the dataframe in case of using soindex or multiple models
)
    inputN = size(f, 1)

    A = Array{Float64}(undef, inputN, inputN)
    for i in 1:inputN
        for j in 1:inputN
            A[i,j] = phi(inputs[i, :],inputs[j, :],k)
        end
    end
    B = [ones(inputN, 1) inputs]
    Bt = transpose(B)


    Mlin = [A B; Bt zeros(size(B, 2),size(B, 2))]

    Flin = [f; zeros(size(B, 2), 1)]

    wv = Mlin \ Flin

    w = wv[1:inputN,:]
    v = wv[inputN+1:size(wv, 1),:]

    polymodel = PolySpline(w, v, inputs, k, orderfunc,name)
    return polymodel
end

function calcpolyspline(
    po::PolySpline,
    x::Array{}
    )
    rv = 0

    for i in 1:size(po.ci, 1)
        fp = po.w[i,1] * phi(x, po.ci[i, :], po.k) # adding radial basis function
        rv += fp[1,1]
    end
    xv = [ones(1,1); x]
    rv += dot(po.v, xv) #adding polynomial term
    return rv
end

function evaluate!(
    po::PolySpline,
    df::DataFrame
    )
    # this function is designed to act like evaluate! in model.jl
    # it takes a model and dataframe and creates a new column: model.name
    # the model is appplied to each row of the dataframe and stashes the solution in the new column
    # use is only possible if the PolySpline has a real orderfunction
    psinput = po.orderfunc(df)
    out = Vector{Float64}(undef, size(psinput, 1))

    for i in 1:size(psinput, 1)
        out[i] = calcpolyspline(po, psinput[i, :])
    end
    df[!, po.name] = out
end
