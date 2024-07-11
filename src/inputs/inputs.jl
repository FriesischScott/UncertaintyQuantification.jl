"""
	sample(inputs::Vector{<:UQInput}, n::Integer)

Generates n correlated samples from a collection of inputs. Returns a DataFrame

See also: [`RandomVariable`](@ref), [`Parameter`](@ref)
"""
function sample(inputs::Vector{<:UQInput}, n::Integer=1)
    return mapreduce(i -> sample(i, n), hcat, inputs)
end

function to_physical_space!(inputs::Vector{<:UQInput}, x::DataFrame)
    for i in inputs
        to_physical_space!(i, x)
    end
    return nothing
end

function to_standard_normal_space!(inputs::Vector{<:UQInput}, x::DataFrame)
    for i in inputs
        to_standard_normal_space!(i, x)
    end
    return nothing
end

function names(inputs::Vector{<:UQInput})
    _names = Symbol[]

    for i in inputs
        if i isa JointDistribution
            append!(_names, names(i))
        else
            push!(_names, getproperty(i, :name))
        end
    end

    return _names
end

function count_rvs(inputs::Vector{<:UQInput})
    random_inputs = filter(i -> (isa(i, RandomUQInput) || isa(i, ProbabilityBox)), inputs)
    return mapreduce(dimensions, +, random_inputs)
end

mean(inputs::Vector{<:UQInput}) = mapreduce(mean, vcat, inputs)

sns_zero_point(inputs::Vector{<:UQInput}) = mapreduce(sns_zero_point, vcat, inputs)
sns_zero_point(input :: UQInput) = (isa(input, RandomUQInput) || isa(input, ProbabilityBox)) ? 0.0 : mean(input)
