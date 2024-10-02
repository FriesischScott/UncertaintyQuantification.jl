"""
	Interval(lb::Real, up::real, name::Symbol)

Defines an Interval, with lower a bound, an upper bound and a name.

# Examples

```jldoctest
julia> Interval(0.10, 0.14, :b)
Interval(0.1, 0.14, :b)```
"""
struct Interval <: ImpreciseUQInput
    lb::Real
    ub::Real
    name::Symbol
    function Interval(lb::Real, ub::Real, name::Symbol)
        lb > ub && error(
            "Lower bound parameter must be smaller than upper bound parameter for Interval $name.",
        )
        return new(lb, ub, name)
    end
end

function map_to_precise(x::Real, input::Interval)
    if !in(x, input)
        error("$x not in [$(input.lb), $(input.ub)] for Interval $(input.name).")
    end
    return Parameter(x, input.name)
end

function sample(i::Interval, n::Integer=1)
    return DataFrame(i.name => fill(i, n))
end

to_standard_normal_space!(i::Interval, x::DataFrame) = nothing
to_physical_space!(i::Interval, x::DataFrame) = nothing

mean(i::Interval) = i

function bounds(i::Interval)
    return i.lb, i.ub
end

Base.in(u, i::Interval) = i.lb <= u <= i.ub

function isdegenerate(i::Interval)
    return i.lb == i.ub
end
