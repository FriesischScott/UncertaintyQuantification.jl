function probability_of_failure(
    models::Union{Vector{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Vector{<:PreciseUQInput},PreciseUQInput},
    sim::AbstractMonteCarlo,
)
    samples = sample(inputs, sim)
    evaluate!(models, samples)

    # Probability of failure
    pf = sum(performance(samples) .< 0) / sim.n

    variance = (pf - pf^2) / sim.n
    cov = sqrt(variance) / pf

    return pf, cov, samples
end

function probability_of_failure(
    models::Union{Vector{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Vector{<:PreciseUQInput},PreciseUQInput},
    sim::LineSampling,
)
    if isempty(sim.direction)
        sim.direction = gradient_in_standard_normal_space(
            [models..., Model(x -> -1 * performance(x), :performance)],
            inputs,
            DataFrame(names(inputs) .=> mean(inputs)),
            :performance,
        )
    end

    samples = sample(inputs, sim)
    evaluate!(models, samples)

    p = reshape(performance(samples), length(sim.points), sim.lines)

    ϕ = Normal()
    ξ = zeros(sim.lines)
    x = median(sim.points)
    for i in 1:(sim.lines)
        if all(p[:, i] .< 0)
            ξ[i] = 1.0
            @warn "All samples for line $i are inside the failure domain"
            continue
        elseif all(p[:, i] .> 0)
            ξ[i] = 0.0
            @warn "All samples for line $i are outside the failure domain"
            continue
        end
        spl = Spline1D(sim.points, p[:, i])
        try
            root = Dierckx.roots(spl)[1]
            ξ[i] = cdf.(ϕ, -root)
        catch e
            @warn "Intersection with failure domain not found for line $i ($e)"
        end
    end

    pf = mean(ξ)
    variance = var(ξ) / sim.lines
    cov = sqrt(variance) / pf

    return pf, cov, samples
end

function probability_of_failure(
    models::Union{Vector{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Vector{<:PreciseUQInput},PreciseUQInput},
    sim::ImportanceSampling,
)
    samples, weights = sample(inputs, sim)
    evaluate!(models, samples)

    # Probability of failure
    weighted_failures = (performance(samples) .< 0) .* weights
    pf = sum(weighted_failures) / sim.n

    variance = ((sum(weighted_failures .* weights) / sim.n) - pf^2) / sim.n
    cov = sqrt(variance) / pf

    return pf, cov, samples
end

# Allow to calculate the pf using only a performance function but no model
function probability_of_failure(
    performance::Function, inputs::Union{Vector{<:PreciseUQInput},PreciseUQInput}, sim::Any
)
    return probability_of_failure(UQModel[], performance, wrap(inputs), sim)
end

# ImpreciseUQInput functions

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

    rhobeg = minimum((ub .- lb) ./ 2)
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

    imprecise_names = [i.name for i in imprecise_inputs]

    g_intervals = map(1:n) do _
        df = sample(precise_inputs)

        bounds = sample.(imprecise_inputs)
        lb = getindex.(bounds, 1)
        ub = getindex.(bounds, 2)

        function g(x)
            inner_sample = hcat(df, DataFrame(imprecise_names .=> x))
            evaluate!(models, inner_sample)

            return performance(inner_sample)[1]
        end

        x0 = (lb .+ ub) ./ 2
        rhobeg = minimum((ub .- lb) ./ 2)

        _, info_min = prima(g, x0; xl=lb, xu=ub, rhobeg=rhobeg)
        g_lb = info_min.fx
        _, info_max = prima(x -> -g(x), x0; xl=lb, xu=ub, rhobeg=rhobeg)
        g_ub = -info_max.fx

        return [g_lb, g_ub]
    end

    pf_ub = sum(getindex.(g_intervals, 1) .< 0) / n
    pf_lb = sum(getindex.(g_intervals, 1) .< 0 .&& getindex.(g_intervals, 2) .< 0) / n

    return Interval(pf_lb, pf_ub, :pf)
end

function bounds(inputs::AbstractVector{UQInput})
    imprecise_inputs = filter(x -> isa(x, ImpreciseUQInput), inputs)

    lb = mapreduce(x -> x.lb, vcat, imprecise_inputs)
    ub = mapreduce(x -> x.ub, vcat, imprecise_inputs)
    return lb, ub
end

function map_to_precise_inputs(x::AbstractVector, inputs::AbstractVector{<:UQInput})
    precise_inputs = PreciseUQInput[]
    params = deepcopy(x)
    for i in inputs
        if isa(i, Interval)
            push!(precise_inputs, map_to_precise(popfirst!(params), i))
        elseif isa(i, ProbabilityBox)
            d = length(i.lb)
            p = [popfirst!(params) for _ in 1:d]
            push!(precise_inputs, map_to_precise(p, i))
        end
    end
    return precise_inputs
end