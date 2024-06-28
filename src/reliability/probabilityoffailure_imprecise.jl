## 1st Method -  External: Global Opt ; Internal: Sampling
function probability_of_failure(
    models::Union{Vector{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
    sim::AbstractSimulation,
)
    inputs = wrap(inputs)
    imprecise_inputs = filter(x -> isa(x, ImpreciseUQInput), inputs)
    precise_inputs = filter(x -> isa(x, PreciseUQInput), inputs)

    function pf(x)
        imprecise_inputs_x = map_to_precise_inputs(x, imprecise_inputs)
        mc_inputs = [precise_inputs..., imprecise_inputs_x...]
        mc_pf, _, _ = probability_of_failure(models, performance, mc_inputs, sim)
        return mc_pf
    end

    lb, ub = float.(bounds(inputs))
    x0 = middle.(lb, ub)

    result_lb = minimize(
        isa(sim, FORM) ? OrthoMADS(length(x0)) : RobustOrthoMADS(length(x0)),
        x -> pf(x),
        x0;
        lowerbound=lb,
        upperbound=ub,
        min_mesh_size=1e-13,
    )

    result_ub = minimize(
        isa(sim, FORM) ? OrthoMADS(length(x0)) : RobustOrthoMADS(length(x0)),
        x -> -pf(x),
        x0;
        lowerbound=lb,
        upperbound=ub,
        min_mesh_size=1e-13,
    )

    pf_lb = result_lb.f
    pf_ub = -result_ub.f

    if pf_lb == pf_ub
        return pf_ub
    else
        return Interval(pf_lb, pf_ub, :pf)
    end
end

## 2nd Method -  External: Sampling ; Internal: Global Opt
function probability_of_failure(
    models::Union{Vector{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
    n::Integer,
)
    inputs = wrap(inputs)
    imprecise_inputs = filter(x -> isa(x, ImpreciseUQInput), inputs)
    precise_inputs = convert(
        Vector{PreciseUQInput}, filter(x -> isa(x, PreciseUQInput), inputs)
    )

    imprecise_names = names(imprecise_inputs)

    g_intervals = map(1:n) do _
        if !isempty(precise_inputs)
            df = sample(precise_inputs)
        else
            df = DataFrame()
        end

        intervals = sample.(imprecise_inputs)
        lb = float.(getindex.(intervals, 1))
        ub = float.(getindex.(intervals, 2))

        function g(x)
            inner_sample = hcat(df, DataFrame(imprecise_names .=> x))
            evaluate!(models, inner_sample)

            return performance(inner_sample)[1]
        end

        x0 = middle.(lb, ub)

        result_lb = minimize(
            OrthoMADS(length(x0)),
            x -> g(x),
            x0;
            lowerbound=lb,
            upperbound=ub,
            min_mesh_size=1e-13,
        )

        result_ub = minimize(
            OrthoMADS(length(x0)),
            x -> -g(x),
            x0;
            lowerbound=lb,
            upperbound=ub,
            min_mesh_size=1e-13,
        )

        return [result_lb.f, -result_ub.f]
    end

    pf_ub = sum(getindex.(g_intervals, 1) .< 0) / n

    pf_lb = sum(getindex.(g_intervals, 2) .< 0) / n

    if pf_lb == pf_ub
        return pf_ub
    else
        return Interval(pf_lb, pf_ub, :pf)
    end
end

function bounds(inputs::AbstractVector{<:UQInput})
    imprecise_inputs = filter(x -> isa(x, ImpreciseUQInput), inputs)

    b = bounds.(imprecise_inputs)
    lb = vcat(getindex.(b, 1)...)
    ub = vcat(getindex.(b, 2)...)
    return lb, ub
end

function map_to_precise_inputs(x::AbstractVector, inputs::AbstractVector{<:UQInput})
    precise_inputs = PreciseUQInput[]
    params = collect(x)
    for i in inputs
        if isa(i, Interval)
            push!(precise_inputs, map_to_precise(popfirst!(params), i))
        elseif isa(i, ProbabilityBox)
            d = length(filter(x -> isa(x, Interval), i.parameters))
            p = [popfirst!(params) for _ in 1:d]
            push!(precise_inputs, map_to_precise(p, i))
        end
    end
    return precise_inputs
end
