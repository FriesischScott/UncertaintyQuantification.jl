struct Model <: AbstractModel
    func::Function
    name::String
end

function (obj::Model)(df::DataFrame)
    obj.func(df)
end

function evaluate(m::Model, df::DataFrame)
    df[!, Symbol(m.name)] = m.func(df)
    return df
end
