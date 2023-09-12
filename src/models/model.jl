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

function evaluate!(models::Vector{<:UQModel}, df::DataFrame)
    # use same working directory for all external models
    datetime = Dates.format(now(), "YYYY-mm-dd-HH-MM-SS")
    for m in models
        if isa(m, ExternalModel)
            evaluate!(m, df; datetime=datetime)
        else
            evaluate!(m, df)
        end
    end
    return nothing
end
