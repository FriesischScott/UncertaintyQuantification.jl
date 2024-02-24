struct Model <: UQModel
    func::Function
    name::Symbol
end

struct ParallelModel <: UQModel
    func::Function
    name::Symbol
end

function (m::Model)(df::DataFrame)
    return m.func(df)
end

function (m::ParallelModel)(df::DataFrame)
    return pmap(m.func, eachrow(df))
end

function evaluate!(m::Model, df::DataFrame)

    if any(eltype.(eachcol(df)) .== Interval)
        return evaluate_imprecise!(m, df)
    end

    return evaluate_precise!(m, df)
end

function evaluate!(m::ParallelModel, df::DataFrame)
    df[!, m.name] = pmap(m.func, eachrow(df))
    return nothing
end

function evaluate_precise!(m::Model, df::DataFrame)
    df[!, m.name] = m.func(df)
    return nothing
end

function evaluate_imprecise!(m::Model, df::DataFrame)

    n = size(df, 1)

    rv_names = names(df)

    imprecise_inputs = eltype.(eachcol(df)) .== Interval
    imprecise_names = rv_names[imprecise_inputs]
    precise_names = rv_names[.~imprecise_inputs]

    g_intervals = fill(Interval(0, 0, m.name), n)

    for i = 1:n

        df_precise = df[i, precise_names]
        df_imprecise = df[i, imprecise_names]

        if isempty(df_precise)
            df_precise = DataFrame()
        else
            df_precise = DataFrame(df_precise)
        end

        bounds = collect(values(df_imprecise))
        lb = getproperty.(bounds, :lb)
        ub = getproperty.(bounds, :ub)

        function f(x)
            df_run = hcat(DataFrame(imprecise_names .=> x), df_precise)

            return m.func(df_run)[1]
        end

        x0 = (lb .+ ub) ./ 2
        rhobeg = minimum(ub .- lb) ./ 4

        _, info_min = prima(f, x0; xl=lb, xu=ub, rhobeg=rhobeg)
        g_lb = info_min.fx
        _, info_max = prima(x -> -f(x), x0; xl=lb, xu=ub, rhobeg=rhobeg)
        g_ub = -info_max.fx

        g_intervals[i] = Interval(g_lb, g_ub, m.name)
    end

    df[!, m.name] .= g_intervals
    return nothing
end

function evaluate!(models::Vector{<:UQModel}, df::DataFrame)
    # use same working directory for all external models
    datetime = Dates.format(now(), "YYYY-mm-dd-HH-MM-SS")
    for m in models
        if isa(m, ExternalModel)
            evaluate!(m, df; datetime=datetime)
        else
            evaluate!(m, df)
        end
    end
    return nothing
end
