mutable struct AdvancedLineSampling <: AbstractSimulation
    lines::Integer
    points::Vector{<:Real}
    direction::NamedTuple
    tolerance::Float64
    stepsize::Float64
    maxiterations::Int64

    function AdvancedLineSampling(
        lines::Integer,
        points::Vector{<:Real}=collect(0.5:0.5:5),
        direction::NamedTuple=NamedTuple(),
        tolerance::Float64=1e-3,
        stepsize::Float64=1e-5,
        maxiterations::Int64=10
    )
        return new(lines, points, direction, tolerance, stepsize, maxiterations)
    end
end

function AdvancedLineSampling(
    lines::Integer,
    tolerance::Float64,
    stepsize::Float64,
    maxiterations::Int64
)
    return AdvancedLineSampling(lines, collect(0.5:0.5:5), NamedTuple(), tolerance, stepsize, maxiterations)
end


function splinefit(
    performance::Function,
    pointsalongline::DataFrame,
    sim::AdvancedLineSampling,
)
    # Evaluate performance function
    p = performance(pointsalongline)

    if all(p[:] .< 0)
        return Inf
    elseif all(p[:] .> 0)
        return -Inf
    end
    spl = Spline1D(sim.points, p[:])
    try
        return Dierckx.roots(spl)[1]
    catch e
        return Inf
    end
end

function _grad(x::Float64, f::Function, sim::AdvancedLineSampling)
    fₓ = only(f(x))
    ∇fₓ = (only(f(x+sim.stepsize)) - fₓ)/sim.stepsize
    return ∇fₓ, fₓ
end


function newtonraphson(x₀::Float64, f::Function, sim::AdvancedLineSampling)
    # Initialize
    err = Inf
    steps = 0

    # Vector to store samples
    x = Float64[]

    while abs(err) > sim.tolerance
        # One dimensional Netwon's method
        ∇fₓ, fₓ = _grad(x₀, f, sim)
        x₁ = x₀ - fₓ/∇fₓ

        append!(x, [x₀, x₀+sim.stepsize])

        # Calculate error
        err = norm(x₁ - x₀)/x₀
        x₀ = x₁

        # Counter
        steps +=1

        if steps >= sim.maxiterations
            @warn "No root found after $steps iterations."
            return Inf
        end
    end

    return x₀, x
end
