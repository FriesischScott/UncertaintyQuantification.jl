function sample(inputs::Array{<:AbstractInput}, n::Int64 = 1)
    mapreduce(i -> sample(i, n), hcat, inputs)
end
