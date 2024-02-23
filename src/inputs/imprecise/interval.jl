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
        lb ≥ ub && error(
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

function sample(i::Interval)
    return [i.lb, i.ub]
end

function bounds(i::Interval)
    return i.lb, i.ub
end

function names(i::AbstractVector{Interval})
    return getproperty.(i, :name)
end

Base.in(u, i::Interval) = i.lb <= u <= i.ub
