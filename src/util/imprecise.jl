function isimprecise(inputs::AbstractVector{<:UQInput})
    return any(isimprecise.(inputs))
end

function isimprecise(input::UQInput)
    return isa(input, Interval) || isa(input, RandomVariable{<:ProbabilityBox})
end
