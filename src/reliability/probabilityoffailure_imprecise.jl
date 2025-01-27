## 1st Method -  External: Global Opt ; Internal: Sampling
struct DoubleLoop
    lb::AbstractSimulation
    ub::AbstractSimulation
end

struct RandomSlicing
    lb::AbstractSimulation
    ub::AbstractSimulation
end

function RandomSlicing(sim::AbstractSimulation)
    return RandomSlicing(sim, deepcopy(sim))
end

function DoubleLoop(sim::AbstractSimulation)
    return DoubleLoop(sim, deepcopy(sim))
end

function probability_of_failure(
    models::Union{Vector{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
    dl::DoubleLoop,
)
    @assert isimprecise(inputs)

    inputs = wrap(inputs)
    imprecise_inputs = filter(x -> isa(x, ImpreciseUQInput), inputs)
    precise_inputs = filter(x -> !isa(x, ImpreciseUQInput), inputs)

    function pf_low(x)
        imprecise_inputs_x = map_to_precise_inputs(x, imprecise_inputs)
        mc_inputs = [precise_inputs..., imprecise_inputs_x...]
        mc_pf, _, _ = probability_of_failure(models, performance, mc_inputs, dl.lb)
        return mc_pf
    end

    function pf_high(x)
        imprecise_inputs_x = map_to_precise_inputs(x, imprecise_inputs)
        mc_inputs = [precise_inputs..., imprecise_inputs_x...]
        mc_pf, _, _ = probability_of_failure(models, performance, mc_inputs, dl.ub)
        return mc_pf
    end

    lb, ub = float.(bounds(inputs))
    x0 = middle.(lb, ub)

    result_lb = minimize(
        isa(dl.lb, FORM) ? OrthoMADS(length(x0)) : RobustOrthoMADS(length(x0)),
        x -> pf_low(x),
        x0;
        lowerbound=lb,
        upperbound=ub,
        min_mesh_size=1e-13,
    )

    result_ub = minimize(
        isa(dl.ub, FORM) ? OrthoMADS(length(x0)) : RobustOrthoMADS(length(x0)),
        x -> -pf_high(x),
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
        return Interval(pf_lb, pf_ub, :pf), result_lb.x, result_ub.x
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

function probability_of_failure(
    models::Union{Vector{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
    rs::RandomSlicing,
)
    @assert isimprecise(inputs)

    inputs = wrap(inputs)

    sns_inputs = mapreduce(transform_to_sns_input, vcat, inputs)

    models = [wrap(models)..., Model(x -> performance(x), :g_slice)]

    sm_min = SlicingModel(models, inputs, false)

    out_ub = probability_of_failure(sm_min, df -> df.g_slice, sns_inputs, rs.ub)

    sm_max = SlicingModel(models, inputs, true)

    out_lb = probability_of_failure(sm_max, df -> df.g_slice, sns_inputs, rs.lb)

    return Interval(out_lb[1], out_ub[1], :pf), out_lb[2:end], out_ub[2:end]
end

function transform_to_sns_input(i::UQInput)
    if isa(i, RandomVariable) || isa(i, ProbabilityBox)
        return RandomVariable(Normal(), i.name)
    elseif isa(i, JointDistribution)
        return RandomVariable.(Normal(), names(i))
    elseif isa(i, Parameter)
        return i
    end

    return UQInput[]
end
