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

function names(inputs::Array{<:UQInput})
    _names = Vector{Symbol}()

    for i in inputs
        if i isa Parameter || i isa RandomVariable
            push!(_names, i.name)
        elseif i isa RandomVariableSet
            append!(_names, names(i))
        end
    end

    return _names
end

mean(inputs::Array{<:UQInput}) = mapreduce(mean, hcat, inputs)
