function isimprecise(inputs::AbstractVector{<:UQInput})
    return any(map(i -> isa(i, ImpreciseUQInput), inputs))
end
