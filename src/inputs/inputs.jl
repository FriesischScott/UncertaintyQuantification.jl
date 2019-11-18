function sample(inputs::Array{<:UQInput}, n::Int64 = 1)
    mapreduce(i -> sample(i, n), hcat, inputs)
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

function _random_inputs(inputs::Array{<:UQInput})
    input_names = Vector{Symbol}()

    for i in inputs
        if isa(i, RandomVariable)
            push!(input_names, i.name)
        elseif isa(i, RandomVariableSet)
            append!(input_names, names(i))
        end
    end

    return input_names
end

mean(inputs::Array{<:UQInput}) = mapreduce(mean, hcat, inputs)
