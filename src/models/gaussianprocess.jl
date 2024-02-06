"""
    GaussianProcess(data::DataFrame, dependendVarName::Symbol, deg::Int, dim::Int)

Creates a gaussian process prior ....

# Examples
```jldoctest
julia> data = DataFrame(x = 1:10, y = [1, 4, 10, 15, 24, 37, 50, 62, 80, 101]);

julia> rs = ResponseSurface(data, :y, 2) |> DisplayAs.withcontext(:compact => true)
ResponseSurface([0.483333, -0.238636, 1.01894], :y, [:x], 2, Monomial{Commutative{CreationOrder}, Graded{LexOrder}}[1, x₁, x₁²])
```
"""
mutable struct GaussianProcessRegressor <: UQModel
    gp::AbstractGP
    y::Symbol
    names::Vector{Symbol}
    function GaussianProcessRegressor(gp::AbstractGP, data::DataFrame, output::Symbol)
        # Choice for kernel and mean function is in gp, is that fine?
        # Where to put output noise?
        names = propertynames(data[:, Not(output)])
        return new(gp, output, names)
    end
end

function fit!(gpr::GaussianProcessRegressor)
    
end



