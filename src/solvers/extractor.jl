struct Extractor
    f::Function
    name::Symbol
end

function names(extractors::Array{Extractor,1})
    return map(e -> e.name, extractors)
end
