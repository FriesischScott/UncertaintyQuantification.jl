struct Model
    func::Function
    name::String
end

function evaluate(m::Model, df::DataFrame)
    m.func(df)
end
