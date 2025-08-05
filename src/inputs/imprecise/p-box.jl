"""
	ProbabilityBox{T}(p::Dict{Symbol,Union{Real,Interval}})

Defines an ProbabilityBox from a `Dict` mapping each of the parameters of the `UnivariateDistribution` `T` to a `Real` or `Interval`.

# Examples

```jldoctest
julia>  ProbabilityBox{Uniform}(Dict(:a => Interval(1.75, 1.83), :b => Interval(1.77, 1.85)))
ProbabilityBox{Uniform}(Dict{Symbol, Union{Real, Interval}}(:a => [1.75, 1.83], :b => [1.77, 1.85]), 1.75, 1.85)
```

```jldoctest
julia>  ProbabilityBox{Normal}(Dict(:μ => Interval(0, 1), :σ =>  Interval(0.1, 1)))
ProbabilityBox{Normal}(Dict{Symbol, Union{Real, Interval}}(:μ => [0, 1], :σ => [0.1, 1]), -Inf, Inf)
```
"""
struct ProbabilityBox{T<:UnivariateDistribution}
    parameters::Dict{Symbol,Union{<:Real,Interval}}
    lb::Real
    ub::Real

    function ProbabilityBox{T}(
        p::Dict{Symbol,Union{Any}}, lb::Real, ub::Real
    ) where {T<:UnivariateDistribution}
        # Make sure all required parameters for the distribution are present
        if !issetequal(keys(p), fieldnames(T))
            error(
                "Parameter mismatch for ProbabilityBox $(keys(p)) != $([fieldnames(T)...])."
            )
        end
        # If someone only passes Parameters, return a RandomVariable instead.
        if all(isa.(values(p), Real))
            @warn "ProbabilityBox() returns a UnivariateDistribution if no intervals are passed"
            return T(getindex.(Ref(p), fieldnames(t))...)
        end
        return new(convert(Dict{Symbol,Union{Real,Interval}}, p), lb, ub)
    end
end

function ProbabilityBox{T}(p::Dict{Symbol,Any}) where {T<:UnivariateDistribution}
    # p-boxes with Uniform distribution as parameter must be treated separately since their support changes with p-box lower and upper bounds.
    if T == Uniform
        parameters = collect(getindex.(Ref(p), fieldnames(T)))
        values = vcat(
            [
                isa(p, Interval) ? collect(UncertaintyQuantification.bounds(p)) : p for
                p in parameters
            ]...,
        )
        return ProbabilityBox{T}(p, minimum(values), maximum(values))
    else
        domain = support(T())
        return ProbabilityBox{T}(p, domain.lb, domain.ub)
    end
end

function ProbabilityBox{T}(parameter::Interval) where {T<:UnivariateDistribution}
    @assert length(fieldnames(T)) == 1
    return ProbabilityBox{T}(Dict{Symbol,Any}(fieldnames(T)[1] => parameter))
end

function ProbabilityBox{T}(p::Dict{Symbol,Interval}) where {T<:UnivariateDistribution}
    return ProbabilityBox{T}(convert(Dict{Symbol,Any}, p))
end

function map_to_precise(
    x::AbstractVector{<:Real}, pbox::ProbabilityBox{T}
) where {T<:UnivariateDistribution}
    parameters = collect(getindex.(Ref(pbox.parameters), fieldnames(T)))
    intervals = filter(x -> isa(x, Interval), parameters)
    if !all(in.(x, intervals))
        error("Values outside of parameter intervals for ProbabilityBox")
    end

    _x = copy(x)

    p = [
        if isa(par, Interval)
            popfirst!(_x)
        else
            par
        end for par in parameters
    ]

    dist_support = support(T())

    if pbox.lb == dist_support.lb && pbox.ub == dist_support.ub
        return T(p...)
    end

    return truncated(T(p...), pbox.lb, pbox.ub)
end

function quantile(pbox::ProbabilityBox{T}, u::Real) where {T<:UnivariateDistribution}
    quantiles = map(
        par -> quantile(map_to_precise([par...], pbox), u),
        Iterators.product([[a, b] for (a, b) in zip(bounds(pbox)...)]...),
    )

    return Interval(minimum(quantiles), maximum(quantiles))
end

rand(pbox::ProbabilityBox, n::Integer=1) = quantile.(Ref(pbox), rand(n))

function bounds(pbox::ProbabilityBox{T}) where {T<:UnivariateDistribution}
    intervals = filter(
        x -> isa(x, Interval), collect(getindex.(Ref(pbox.parameters), fieldnames(T)))
    )
    lb = getproperty.(intervals, :lb)
    ub = getproperty.(intervals, :ub)

    return lb, ub
end

function cdf(pbox::ProbabilityBox{T}, x::Real) where {T<:UnivariateDistribution}
    lb, ub = bounds(pbox)

    cdfs_lo = map(
        par -> cdf(map_to_precise([par...], pbox), x),
        Iterators.product([[a, b] for (a, b) in zip(lb, ub)]...),
    )

    cdfs_hi = map(
        par -> cdf(map_to_precise([par...], pbox), x),
        Iterators.product([[a, b] for (a, b) in zip(lb, ub)]...),
    )

    return Interval(minimum(cdfs_lo), maximum(cdfs_hi))
end

# Does the inverse of quantile, not cdf, which would return an interval
function reverse_quantile(
    pbox::ProbabilityBox{T}, x::Interval
) where {T<:UnivariateDistribution}
    lb, ub = bounds(pbox)

    cdfs_lo = map(
        par -> cdf(map_to_precise([par...], pbox), x.lb),
        Iterators.product([[a, b] for (a, b) in zip(lb, ub)]...),
    )

    cdfs_hi = map(
        par -> cdf(map_to_precise([par...], pbox), x.ub),
        Iterators.product([[a, b] for (a, b) in zip(lb, ub)]...),
    )

    u_lo = maximum(cdfs_lo)
    u_hi = minimum(cdfs_hi)

    error = abs(u_hi - u_lo)
    if error > 10^-6
        @warn(
            "When inverting the quantile function for a p-box, the error was $(error), greater than the allowed tolerance of 10^-6"
        )
    end
    return middle(u_lo, u_hi)   # Return midpoint
end

Base.broadcastable(pbox::ProbabilityBox) = Ref(pbox)

length(::ProbabilityBox{T}) where {T<:UnivariateDistribution} = 1
