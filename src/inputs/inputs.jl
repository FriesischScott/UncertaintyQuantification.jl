function sample(inputs::Array{<:AbstractInput}, n::Int64 = 1)
    mapreduce(i -> sample(i, n), hcat, inputs)
end

function to_physical_space!(inputs::Array{<:AbstractInput}, x::DataFrame)
    for i in inputs
        to_physical_space!(i, x)
    end
    return nothing
end

function to_standard_normal_space!(inputs::Array{<:AbstractInput}, x::DataFrame)
    for i in inputs
        to_standard_normal_space!(i, x)
    end
    return nothing
end

function _random_inputs(inputs::Array{<:AbstractInput})
    names = Vector{Symbol}()

    for i in inputs
        if isa(i, RandomVariable)
            push!(names, Symbol(i.name))
        elseif isa(i, RandomVariableSet)
            append!(names, map(x -> Symbol(x.name), i.members))
        end
    end

    return names
end

mean(inputs::Array{<:AbstractInput}) = mapreduce(mean, hcat, inputs)
