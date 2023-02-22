"""
	sample(inputs::Array{<:UQInput}, n::Integer)

Generates n correlated samples from a collection of inputs. Returns a DataFrame

See also: [`RandomVariable`](@ref), [`Parameter`](@ref)
"""
function sample(inputs::Array{<:UQInput}, n::Integer=1)
    return mapreduce(i -> sample(i, n), hcat, inputs)
end

function to_physical_space!(inputs::Array{<:UQInput}, x::DataFrame)
    for i in inputs
        to_physical_space!(i, x)
    end
    return nothing
end

function to_standard_normal_space!(inputs::Array{<:UQInput}, x::DataFrame)
    for i in inputs
        to_standard_normal_space!(i, x)
    end
    return nothing
end

function names(inputs::Array{<:UQInput})
    _names = Symbol[]

    for i in inputs
        if i isa Parameter || i isa RandomVariable
            push!(_names, i.name)
        elseif i isa JointDistribution
            append!(_names, names(i))
        end
    end

    return _names
end

function count_rvs(inputs::Array{<:UQInput})
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    return mapreduce(dimensions, +, random_inputs)
end

mean(inputs::Array{<:UQInput}) = mapreduce(mean, vcat, inputs)
