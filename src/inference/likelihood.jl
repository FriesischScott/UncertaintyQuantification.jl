struct LogLikelihood <: AbstractLikelihood
    model::UQModel
    data::DataFrame
    distance::Function

    LogLikelihood(model,distance) = new(model,distance)
end

