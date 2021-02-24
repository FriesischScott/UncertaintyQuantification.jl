struct Extractor
    f::Function
    name::Symbol
end

function names(extractors::Array{Extractor,1})
    map(e -> e.name, extractors)
end