# TODO: Add interval propagation method
struct SlicingModel <: UQModel
    models::AbstractVector{<:UQModel}
    inputs::AbstractVector{<:UQInput}
    max::Bool
end

"""
 evaluates the transformed SNS problem
"""
function evaluate!(m::SlicingModel, df::DataFrame)
    pboxes = filter(x -> isa(x, ProbabilityBox), m.inputs)
    intervals = filter(x -> isa(x, Interval), m.inputs)

    if isempty(intervals)
        intervals = Interval[]
    end

    pbox_names = names(pboxes)
    interval_names = names(intervals)

    imprecise_names = [pbox_names..., interval_names...]

    perf = m.models[end].name

    physical = copy(df)

    to_physical_space!(m.inputs, physical)

    g = map(eachrow(physical)) do row
        lb = [getindex.(collect(row[pbox_names]), 1)..., getproperty.(intervals, :lb)...]

        ub = [getindex.(collect(row[pbox_names]), 2)..., getproperty.(intervals, :ub)...]

        x0 = middle.(lb, ub)

        precise_df = hcat(
            DataFrame(
                [[zero(eltype(lb))] for _ in 1:length(imprecise_names)], imprecise_names
            ),
            select(DataFrame(row), Not(imprecise_names)),
        )

        function f(x)
            precise_df[1, imprecise_names] .= x

            evaluate!(m.models, precise_df)

            return precise_df[1, perf]
        end

        result = minimize(
            OrthoMADS(length(x0)),
            x -> m.max ? -f(x) : f(x),
            x0;
            lowerbound=lb,
            upperbound=ub,
            min_mesh_size=1e-13,
        )

        return m.max ? -result.f : result.f
    end

    df[!, perf] = g

    return nothing
end
