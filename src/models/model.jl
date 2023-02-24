struct Model <: UQModel
    func::Function
    name::Symbol
end

struct ParallelModel <: UQModel
    func::Function
    name::Symbol
end

function (m::Model)(df::DataFrame)
    return m.func(df)
end

function (m::ParallelModel)(df::DataFrame)
    return pmap(m.func, eachrow(df))
end

function evaluate!(m::Model, df::DataFrame)
    df[!, m.name] = m.func(df)
    return nothing
end

function evaluate!(m::ParallelModel, df::DataFrame)
    df[!, m.name] = pmap(m.func, eachrow(df))
    return nothing
end

function evaluate!(models::Vector{T} where {T<:UQModel}, df::DataFrame)
    for m in models
        evaluate!(m, df)
    end
    return nothing
end
