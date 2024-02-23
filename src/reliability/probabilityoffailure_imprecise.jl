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

    function montecarlo_pf(x)
        imprecise_inputs_x = map_to_precise_inputs(x, imprecise_inputs)
        mc_inputs = [precise_inputs..., imprecise_inputs_x...]
        mc_pf, _, _ = probability_of_failure(models, performance, mc_inputs, sim)
        return mc_pf
    end

    lb, ub = bounds(inputs)
    x0 = (lb .+ ub) ./ 2

    rhobeg = minimum(ub .- lb) ./ 4
    _, info_min = prima(montecarlo_pf, x0; xl=lb, xu=ub, rhobeg=rhobeg)
    pf_lb = info_min.fx
    _, info_max = prima(x -> -montecarlo_pf(x), x0; xl=lb, xu=ub, rhobeg=rhobeg)
    pf_ub = -info_max.fx

    return Interval(pf_lb, pf_ub, :pf)
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

        bounds = sample.(imprecise_inputs)
        lb = getindex.(bounds, 1)
        ub = getindex.(bounds, 2)

        function g(x)
            inner_sample = hcat(df, DataFrame(imprecise_names .=> x))
            evaluate!(models, inner_sample)

            return performance(inner_sample)[1]
        end

        x0 = (lb .+ ub) ./ 2
        rhobeg = minimum(ub .- lb) ./ 4

        _, info_min = prima(g, x0; xl=lb, xu=ub, rhobeg=rhobeg)
        g_lb = info_min.fx
        _, info_max = prima(x -> -g(x), x0; xl=lb, xu=ub, rhobeg=rhobeg)
        g_ub = -info_max.fx

        return [g_lb, g_ub]
    end

    left_bound_failures = getindex.(g_intervals, 1) .< 0
    pf_ub = sum(left_bound_failures) / n

    both_bound_failures = getindex.(g_intervals[left_bound_failures], 2) .< 0
    pf_lb = sum(both_bound_failures) / n

    return Interval(pf_lb, pf_ub, :pf)
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
    params = copy(x)
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
