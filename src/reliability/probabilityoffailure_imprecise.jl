## 1st Method -  External: Global Opt ; Internal: Sampling
struct DoubleLoop
    sim::AbstractSimulation
end

struct IntervalMonteCarlo
    sim::AbstractSimulation
end

function probability_of_failure(
    models::Union{Vector{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
    dl::DoubleLoop,
)
    inputs = wrap(inputs)
    imprecise_inputs = filter(x -> isa(x, ImpreciseUQInput), inputs)
    precise_inputs = filter(x -> !isa(x, ImpreciseUQInput), inputs)

    function pf(x)
        imprecise_inputs_x = map_to_precise_inputs(x, imprecise_inputs)
        mc_inputs = [precise_inputs..., imprecise_inputs_x...]
        mc_pf, _, _ = probability_of_failure(models, performance, mc_inputs, dl.sim)
        return mc_pf
    end

    lb, ub = float.(bounds(inputs))
    x0 = middle.(lb, ub)

    result_lb = minimize(
        isa(dl.sim, FORM) ? OrthoMADS(length(x0)) : RobustOrthoMADS(length(x0)),
        x -> pf(x),
        x0;
        lowerbound=lb,
        upperbound=ub,
        min_mesh_size=1e-13,
    )

    result_ub = minimize(
        isa(dl.sim, FORM) ? OrthoMADS(length(x0)) : RobustOrthoMADS(length(x0)),
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

hi(x::Interval) = x.ub
lo(x::Interval) = x.lb
hi(x::Float64) = x
lo(x::Float64) = x

## 2nd Method -  External: Sampling ; Internal: Global Opt
function probability_of_failure_slicing(
    models::Union{Vector{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
    imc::IntervalMonteCarlo,
)

    model_imprecise = Model(x -> performance(hcat(DataFrame(models.name => models.func(x)), x)), :performance)
    performance_hi(df) = hi.(df[!, :performance])
    performance_lo(df) = lo.(df[!, :performance])

    pf_ub, sig_ub, samples_hi = probability_of_failure(model_imprecise, performance_lo, inputs, imc.sim)
    pf_lb, sig_lb, samples_lo = probability_of_failure(model_imprecise, performance_hi, inputs, imc.sim)

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
    precise_inputs = UQInput[]
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
