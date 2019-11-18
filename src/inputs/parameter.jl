struct Parameter <: DeterministicUQInput
    value::Real
    name::Symbol
end

function sample(p::Parameter, n::Int64 = 1)
    DataFrame(p.name => ones(n) * p.value)
end

to_standard_normal_space!(p::Parameter, x::DataFrame) = nothing
to_physical_space!(p::Parameter, x::DataFrame) = nothing

mean(p::Parameter) = DataFrame(p.name => p.value)
