struct Parameter <: AbstractInput
    value::Real
    name::String
end

function sample(p::Parameter, n::Int64 = 1)
    DataFrame(Symbol(p.name) => ones(n) * p.value)
end
