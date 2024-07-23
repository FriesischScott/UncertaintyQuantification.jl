struct StochasticProcessModel <: UQModel
    proc::AbstractStochasticProcess
    name::Symbol
end

function StochasticProcessModel(proc::AbstractStochasticProcess)
    return StochasticProcessModel(proc, proc.name)
end

function evaluate!(m::StochasticProcessModel, df::DataFrame)
    timeseries = map(eachrow(df)) do row
        return m.proc(row)
    end

    return df[!, m.name] .= timeseries
end
