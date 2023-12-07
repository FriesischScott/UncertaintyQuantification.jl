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
end

function map_to_precise(x::Real, input::Interval)
    x < input.lb && error("Choosen value $x is lower than Interval's lower bound $input.lb")
    x > input.ub && error("Choosen value $x is higher than Interval's upper bound $input.lb")
    return Parameter(x, input.name)
end