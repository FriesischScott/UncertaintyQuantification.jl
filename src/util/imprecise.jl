function isimprecise(inputs::AbstractVector{<:UQInput})
    return any(isimprecise.(inputs))
end

function isimprecise(input::UQInput)
    return isa(input, IntervalVariable) ||
           isa(input, RandomVariable{<:ProbabilityBox}) ||
           (
               isa(input, JointDistribution) &&
               any(isa.(input.marginals, RandomVariable{<:ProbabilityBox}))
           )
end
