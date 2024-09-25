"""
	ProbabilityBox{T}(p::AbstractVector{Interval}, name::Symbol)

Defines an ProbabilityBox from a `Vector` of `Interval`, name `UnivariateDistribution` `T`. The number and order of parameters must match the parameters of the associated distribution from Distributions.jl.

# Examples

```jldoctest
julia>  ProbabilityBox{Uniform}([Interval(1.75, 1.83, :a), Interval(1.77, 1.85, :b)], :l)
ProbabilityBox{Uniform}(Interval[Interval(1.75, 1.83, :a), Interval(1.77, 1.85, :b)], :l, 0.0, 1.0)
```

```jldoctest
julia>  ProbabilityBox{Normal}([Interval(0, 1, :μ), Interval(0.1, 1, :σ)], :x)
ProbabilityBox{Normal}(Interval[Interval(0, 1, :μ), Interval(0.1, 1, :σ)], :x, -Inf, Inf)
```
"""
struct ProbabilityBox{T<:UnivariateDistribution} <: ImpreciseUQInput
    parameters::AbstractVector{<:UQInput}
    name::Symbol
    lb::Real
    ub::Real

    function ProbabilityBox{T}(
        p::AbstractVector{<:UQInput}, name::Symbol, lb::Real, ub::Real
    ) where {T<:UnivariateDistribution}
        # Only allow Intervals and Parameters
        if !isempty(filter(x -> !isa(x, Interval) && !isa(x, Parameter), p))
            error("A ProbabilityBox can only be constructed from Intervals and Parameters.")
        end

        # Make sure all required parameters for the distribution are present
        if !(names(p) == [fieldnames(T)...])
            error(
                "Parameter mismatch for ProbabilityBox $name: $(names(p)) != $([fieldnames(T)...]).",
            )
        end

        # If someone only passes Parameters, return a RandomVariable instead.
        if all(isa.(p, Parameter))
            return RandomVariable(T(getproperty.(p, :value)...), name)
        end

        return new(p, name, lb, ub)
    end
end

function ProbabilityBox{T}(
    p::AbstractVector{<:UQInput}, name::Symbol
) where {T<:UnivariateDistribution}
    if isa(T(), Uniform)
        v_ip = mapreduce(
            x -> collect(bounds(x)), vcat, filter(x -> isa(x, ImpreciseUQInput), p)
        )
        v_p = map(x -> x.value, filter(x -> isa(x, Parameter), p))
        values = vcat(v_ip, v_p)
        return ProbabilityBox{T}(p, name, minimum(values), maximum(values))
    else
        domain = support(T())
        return ProbabilityBox{T}(p, name, domain.lb, domain.ub)
    end
end

function ProbabilityBox{T}(
    parameter::Interval, name::Symbol
) where {T<:UnivariateDistribution}
    return ProbabilityBox{T}([parameter], name)
end

function map_to_precise(
    x::AbstractVector{<:Real}, pbox::ProbabilityBox{T}
) where {T<:UnivariateDistribution}
    intervals = filter(x -> isa(x, Interval), pbox.parameters)
    if !all(in.(x, intervals))
        error("Values outside of parameter intervals for ProbabilityBox $(pbox.name).")
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
        return RandomVariable(T(p...), pbox.name)
    end

    return RandomVariable(truncated(T(p...), pbox.lb, pbox.ub), pbox.name)
end

function quantile(pbox::ProbabilityBox{T}, u::Real) where {T<:UnivariateDistribution}
    lb, ub = bounds(pbox)

    quantiles = map(
        par -> quantile(map_to_precise([par...], pbox), u),
        Iterators.product([[a, b] for (a, b) in zip(lb, ub)]...),
    )

    lb = minimum(quantiles)
    ub = maximum(quantiles)

    if lb == ub
        return lb
    else
        return Interval(lb, ub, pbox.name)
    end
end

rand(pbox::ProbabilityBox, n::Integer=1) = quantile.(Ref(pbox), rand(n))

function bounds(pbox::ProbabilityBox{T}) where {T<:UnivariateDistribution}
    intervals = filter(x -> isa(x, Interval), pbox.parameters)
    lb = getproperty.(intervals, :lb)
    ub = getproperty.(intervals, :ub)

    return lb, ub
end

function sample(pbox::ProbabilityBox, n::Integer=1)
    return DataFrame(pbox.name => rand(pbox, n))
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

    return Interval(minimum(cdfs_lo), maximum(cdfs_hi), :cdf)
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
            "When inverting the quantile function for p-box $(pbox.name), the error was $(error), greater than the allowed tolerance of 10^-6"
        )
    end
    return mean([u_lo, u_hi])   # Return midpoint
end

function to_physical_space!(
    pbox::ProbabilityBox{T}, x::DataFrame
) where {T<:UnivariateDistribution}
    x[!, pbox.name] = map(x -> quantile(pbox, x), cdf.(Normal(), collect(x[:, pbox.name])))
    return nothing
end

function to_standard_normal_space!(
    pbox::ProbabilityBox{T}, x::DataFrame
) where {T<:UnivariateDistribution}
    x[!, pbox.name] = _to_standard_normal_space(pbox, x[:, pbox.name])
    return nothing
end

function _to_standard_normal_space(
    d::ProbabilityBox{T}, x::Vector
) where {T<:UnivariateDistribution}
    return quantile.(Normal(), reverse_quantile.(Ref(d), x))
end

dimensions(d::ProbabilityBox{T}) where {T<:UnivariateDistribution} = 1

length(::ProbabilityBox{T}) where {T<:UnivariateDistribution} = 1
