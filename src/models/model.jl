struct Model <: UQModel
    func::Function
    name::Symbol
end

function (obj::Model)(df::DataFrame)
    obj.func(df)
end

function evaluate(m::Model, df::DataFrame)
    df[!, m.name] = m.func(df)
    return df
end
