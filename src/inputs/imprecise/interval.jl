"""
	Interval(lb::Real, up::real)

Defines an Interval, with lower a bound, an upper bound.

This is a data type used internally and to construct p-boxes. For interval inputs see [`IntervalVariable`](@ref).

# Examples

```jldoctest
julia> Interval(0.10, 0.14)
Interval(0.1, 0.14)```
"""
struct Interval
    lb::Real
    ub::Real
    function Interval(lb::Real, ub::Real)
        lb > ub && error(
            "Lower bound parameter must be smaller than upper bound parameter for Intervals.",
        )
        return new(lb, ub)
    end
end

function Base.show(io::IO, i::Interval)
    print(io, "[$(i.lb), $(i.ub)]")
    return nothing
end

Base.in(u, i::Interval) = i.lb <= u <= i.ub

Base.broadcastable(i::Interval) = Ref(i)

function isdegenerate(i::Interval)
    return i.lb == i.ub
end

function bounds(i::Interval)
    return i.lb, i.ub
end

"""
	IntervalVariable(lb::Real, up::real, name::Symbol)

Defines an IntervalVariable input , with lower bound `lb`, upper bound `ub` and `name`.

This is an input type. To construct p-boxes with interval parameters use [`Interval`](@ref)

# Examples

```jldoctest
julia> Interval(0.10, 0.14, :x)
Interval(0.1, 0.14, :x)```
"""
struct IntervalVariable <: UQInput
    lb::Real
    ub::Real
    name::Symbol
    function IntervalVariable(lb::Real, ub::Real, name::Symbol)
        lb > ub && error(
            "Lower bound parameter must be smaller than upper bound parameter for Interval $name.",
        )
        return new(lb, ub, name)
    end
end

function Base.show(io::IO, i::IntervalVariable)
    print(io, "$(String(i.name)) ∈ [$(i.lb),$(i.ub)]")
    return nothing
end

function Interval(i::IntervalVariable)
    return Interval(i.lb, i.ub)
end

function sample(i::IntervalVariable, n::Integer=1)
    return DataFrame(i.name => fill(Interval(i), n))
end

function map_to_precise(x::Real, input::IntervalVariable)
    if !in(x, input)
        error("$x not in [$(input.lb), $(input.ub)] for Interval $(input.name).")
    end
    return Parameter(x, input.name)
end

to_standard_normal_space!(_::IntervalVariable, _::DataFrame) = nothing
to_physical_space!(_::IntervalVariable, _::DataFrame) = nothing

function bounds(i::IntervalVariable)
    return i.lb, i.ub
end

Base.in(u, i::IntervalVariable) = i.lb <= u <= i.ub

mean(i::IntervalVariable) = Interval(i)
