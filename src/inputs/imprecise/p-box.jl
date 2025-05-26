"""
	ProbabilityBox{T}(p::AbstractVector{<:UQInput})

Defines an ProbabilityBox from a `Vector` of `UQInput` and `UnivariateDistribution` `T`. The number and order of parameters must match the parameters of the associated distribution from Distributions.jl.

# Examples

```jldoctest
julia>  ProbabilityBox{Uniform}([Interval(1.75, 1.83, :a), Interval(1.77, 1.85, :b)], :l)
ProbabilityBox{Uniform}(Interval[Interval(1.75, 1.83, :a), Interval(1.77, 1.85, :b)], :l, 1.75, 1.85)
```

```jldoctest
julia>  ProbabilityBox{Normal}([Interval(0, 1, :μ), Interval(0.1, 1, :σ)], :x)
ProbabilityBox{Normal}(Interval[Interval(0, 1, :μ), Interval(0.1, 1, :σ)], :x, -Inf, Inf)
```
"""
struct ProbabilityBox{T<:UnivariateDistribution}
    parameters::AbstractVector{<:UQInput}
    lb::Real
    ub::Real

    function ProbabilityBox{T}(
        p::AbstractVector{<:UQInput}, lb::Real, ub::Real
    ) where {T<:UnivariateDistribution}
        # Only allow Intervals and Parameters
        if !isempty(filter(x -> !isa(x, Interval) && !isa(x, Parameter), p))
            error("A ProbabilityBox can only be constructed from Intervals and Parameters.")
        end
        # Make sure all required parameters for the distribution are present
        if !(names(p) == [fieldnames(T)...])
            error(
                "Parameter mismatch for ProbabilityBox $(names(p)) != $([fieldnames(T)...]).",
            )
        end
        # If someone only passes Parameters, return a RandomVariable instead.
        if all(isa.(p, Parameter))
            @warn "ProbabilityBox() returns UnivariateDistribution if only Parameters are passed"
            return T(getproperty.(p, :value)...)
        end
        return new(p, lb, ub)
    end
end

function ProbabilityBox{T}(p::AbstractVector{<:UQInput}) where {T<:UnivariateDistribution}
    # p-boxes with Uniform distribution as parameter must be treated separately since their support changes with p-box lower and upper bounds.
    if T == Uniform
        bounds_intervals = mapreduce(x -> collect(bounds(x)), vcat, p[isimprecise.(p)])  # collecting bounds of the intervals used for describing the Uniform distribution
        values_parameters = map(x -> x.value, p[.!isimprecise.(p)]) # collecting values of Parameters used for describing the Uniform distribution
        values = vcat(bounds_intervals, values_parameters)
        return ProbabilityBox{T}(p, minimum(values), maximum(values))
    else
        domain = support(T())
        return ProbabilityBox{T}(p, domain.lb, domain.ub)
    end
end

function ProbabilityBox{T}(parameter::Interval) where {T<:UnivariateDistribution}
    return ProbabilityBox{T}([parameter])
end

function map_to_distribution(
    x::AbstractVector{<:Real}, pbox::ProbabilityBox{T}
) where {T<:UnivariateDistribution}
    intervals = filter(x -> isa(x, Interval), pbox.parameters)
    if !all(in.(x, intervals))
        error("Values outside of parameter intervals for ProbabilityBox")
    end

    _x = copy(x)

    p = [
        if isa(par, Interval)
            popfirst!(_x)
        else
            par.value
        end for par in pbox.parameters
    ]

    dist_support = support(T())

    if pbox.lb == dist_support.lb && pbox.ub == dist_support.ub
        return T(p...)
    end

    return truncated(T(p...), pbox.lb, pbox.ub)
end

function quantile(pbox::ProbabilityBox{T}, u::Real) where {T<:UnivariateDistribution}
    quantiles = map(
        par -> quantile(map_to_distribution([par...], pbox), u),
        Iterators.product([[a, b] for (a, b) in zip(bounds(pbox)...)]...),
    )

    return (lb=minimum(quantiles), ub=maximum(quantiles))
end

rand(pbox::ProbabilityBox, n::Integer=1) = quantile.(Ref(pbox), rand(n))

function bounds(pbox::ProbabilityBox{T}) where {T<:UnivariateDistribution}
    intervals = filter(x -> isa(x, Interval), pbox.parameters)
    lb = getproperty.(intervals, :lb)
    ub = getproperty.(intervals, :ub)

    return lb, ub
end

function cdf(pbox::ProbabilityBox{T}, x::Real) where {T<:UnivariateDistribution}
    lb, ub = bounds(pbox)

    cdfs_lo = map(
        par -> cdf(map_to_distribution([par...], pbox), x),
        Iterators.product([[a, b] for (a, b) in zip(lb, ub)]...),
    )

    cdfs_hi = map(
        par -> cdf(map_to_distribution([par...], pbox), x),
        Iterators.product([[a, b] for (a, b) in zip(lb, ub)]...),
    )

    return Interval(minimum(cdfs_lo), maximum(cdfs_hi), :cdf)
end

# Does the inverse of quantile, not cdf, which would return an interval
function reverse_quantile(
    pbox::ProbabilityBox{T}, x::NamedTuple
) where {T<:UnivariateDistribution}
    lb, ub = bounds(pbox)

    cdfs_lo = map(
        par -> cdf(map_to_distribution([par...], pbox), x.lb),
        Iterators.product([[a, b] for (a, b) in zip(lb, ub)]...),
    )

    cdfs_hi = map(
        par -> cdf(map_to_distribution([par...], pbox), x.ub),
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
