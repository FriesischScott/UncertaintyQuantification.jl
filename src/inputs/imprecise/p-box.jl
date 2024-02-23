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

function quantile(pbox::ProbabilityBox{T}, u::Real) where {T<:UnivariateDistribution}
    lb = getproperty.(pbox.parameters, :lb)
    ub = getproperty.(pbox.parameters, :ub)

    quantiles = map(
        par -> quantile(T(par...), u),
        Iterators.product([[a, b] for (a, b) in zip(lb, ub)]...),
    )

    lb = minimum(quantiles)
    ub = maximum(quantiles)

    return Interval(lb, ub, pbox.name)
end

rand(pbox::ProbabilityBox, n::Integer=1) = quantile.(Ref(pbox), rand(n))

function bounds(pbox::ProbabilityBox{T}) where {T<:UnivariateDistribution}
    lb = getproperty.(pbox.parameters, :lb)
    ub = getproperty.(pbox.parameters, :ub)

    return lb, ub
end

function sample(pbox::ProbabilityBox, n::Integer=1)
    return DataFrame(pbox.name => rand(pbox, n))
end

# Does the inverse of quantile, not cdf, which would return an interval
function reverse_quantile(pbox::ProbabilityBox{T}, x::Interval) where {T<:UnivariateDistribution}
    lb = getproperty.(pbox.parameters, :lb)
    ub = getproperty.(pbox.parameters, :ub)

    cdfs_lo = map(
        par -> cdf(T(par...), x.lb),
        Iterators.product([[a, b] for (a, b) in zip(lb, ub)]...),
    )

    cdfs_hi = map(
        par -> cdf(T(par...), x.ub),
        Iterators.product([[a, b] for (a, b) in zip(lb, ub)]...),
    )

    u_lo = maximum(cdfs_lo)
    u_hi = minimum(cdfs_hi)
    
    error = abs(u_hi - u_lo)
    if error >10^-6
        @warn("When inverting the quantile function for p-box $(pbox.name), the error was $(error)")
    end
    return mean([u_lo, u_hi])   # Return midpoint
end

function to_physical_space!(pbox::ProbabilityBox{T}, x::DataFrame)  where {T<:UnivariateDistribution}
    x[!, pbox.name] = _to_physical_space(pbox, x[:, pbox.name])
    return nothing
end

function _to_physical_space(d::ProbabilityBox{T}, x::Vector) where {T<:UnivariateDistribution}
    return quantile.(Ref(d), cdf.(Normal(), x))
end

function to_standard_normal_space!(pbox::ProbabilityBox{T}, x::DataFrame) where {T<:UnivariateDistribution}
    x[!, pbox.name] = _to_standard_normal_space(pbox, x[:, pbox.name])
    return nothing
end

function _to_standard_normal_space(d::ProbabilityBox{T}, x::Vector) where {T<:UnivariateDistribution}
    return quantile.(Normal(), reverse_quantile.(Ref(d), x))
end

dimensions(d::ProbabilityBox{T}) where {T<:UnivariateDistribution} = 1

length(::ProbabilityBox{T}) where {T<:UnivariateDistribution} = 1