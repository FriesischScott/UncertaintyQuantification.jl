"""
	ProbabilityBox{T}(p::AbstractVector{Interval}, name::Symbol)

Defines an ProbabilityBox from a `Vector` of `Interval`, name `UnivariateDistribution` `T`. The number and order of parameters must match the parameters of the associated distribution from Distributions.jl.

# Examples

```jldoctest
julia>  ProbabilityBox{Uniform}([Interval(1.75, 1.83, :a), Interval(1.77, 1.85, :b)], :l)
ProbabilityBox{Uniform}(Interval[Interval(1.75, 1.83, :a), Interval(1.77, 1.85, :b)], :l)
```

```jldoctest
julia>  ProbabilityBox{Normal}([Interval(0, 1, :μ), Interval(0.1, 1, :σ)], :x)
ProbabilityBox{Normal}(Interval[Interval(0, 1, :μ), Interval(0.1, 1, :σ)], :x)
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
    domain = support(T())
    return ProbabilityBox{T}(p, name, domain.lb, domain.ub)
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

function sample(pbox::ProbabilityBox{T}, u::Real) where {T<:UnivariateDistribution}
    lb, ub = bounds(pbox)

    quantiles = map(
        par -> quantile(map_to_precise([par...], pbox), u),
        Iterators.product([[a, b] for (a, b) in zip(lb, ub)]...),
    )

    return [minimum(quantiles), maximum(quantiles)]
end

sample(pbox::ProbabilityBox{T}) where {T<:UnivariateDistribution} = sample(pbox, rand())

function bounds(pbox::ProbabilityBox{T}) where {T<:UnivariateDistribution}
    intervals = filter(x -> isa(x, Interval), pbox.parameters)
    lb = getproperty.(intervals, :lb)
    ub = getproperty.(intervals, :ub)

    return lb, ub
end
