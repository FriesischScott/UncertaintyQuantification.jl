struct Extractor
    f::Function
    name::Symbol
end

function names(extractors::Vector{Extractor})
    return map(e -> e.name, extractors)
end
