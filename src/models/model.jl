struct Model <: UQModel
    func::Function
    name::Symbol
end

function (obj::Model)(df::DataFrame)
    return obj.func(df)
end

function evaluate!(models::Array{T,1} where {T<:UQModel}, df::DataFrame)
    for m in models
        evaluate!(m, df)
    end
    return nothing
end

function evaluate!(m::Model, df::DataFrame)
    df[!, m.name] = m.func(df)
    return nothing
end
