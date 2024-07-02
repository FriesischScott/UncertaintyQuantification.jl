wrap(x) = [x]
wrap(x::AbstractArray) = x
wrap(x::Vector{<:UQInput}) = x
