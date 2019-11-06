function sample(inputs::Array{<:AbstractInput}, n::Int64 = 1)
    map(i -> sample(i, n), inputs) |> samples -> hcat(samples...)
end
