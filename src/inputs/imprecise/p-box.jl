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
ProbabilityBox{Normal}([Interval(0, 1, :μ), Interval(0.1, 1, :σ)], :x)
```
"""
struct ProbabilityBox{T<:UnivariateDistribution} <: ImpreciseUQInput
    parameters::AbstractVector{Interval}
    name::Symbol

    function ProbabilityBox{T}(
        p::AbstractVector{Interval}, name::Symbol
    ) where {T<:UnivariateDistribution}
        if !(names(p) == [fieldnames(T)...])
            error(
                "Parameter mismatch for ProbabilityBox $name: $(names(p)) != $([fieldnames(T)...]).",
            )
        end

        return new(p, name)
    end
end

function ProbabilityBox{T}(
    parameter::Interval, name::Symbol
) where {T<:UnivariateDistribution}
    return ProbabilityBox{T}([parameter], name)
end

function map_to_precise(
    x::Vector{<:Real}, pbox::ProbabilityBox{T}
) where {T<:UnivariateDistribution}
    if !all(in.(x, pbox.parameters))
        error("Values outside of parameter intervals for ProbabilityBox $(pbox.name).")
    end

    return RandomVariable(T(x...), pbox.name)
end

function sample(pbox::ProbabilityBox{T}, u::Real) where {T<:UnivariateDistribution}
    lb = getproperty.(pbox.parameters, :lb)
    ub = getproperty.(pbox.parameters, :ub)

    quantiles = map(
        par -> quantile(T(par...), u),
        Iterators.product([[a, b] for (a, b) in zip(lb, ub)]...),
    )

    return [minimum(quantiles), maximum(quantiles)]
end

sample(pbox::ProbabilityBox{T}) where {T<:UnivariateDistribution} = sample(pbox, rand())

function bounds(pbox::ProbabilityBox{T}) where {T<:UnivariateDistribution}
    lb = getproperty.(pbox.parameters, :lb)
    ub = getproperty.(pbox.parameters, :ub)

    return lb, ub
end
