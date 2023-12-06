using UncertaintyQuantification
using PRIMA

struct Interval <: ImpreciseUQInput
    lb::Real
    ub::Real
    name::Symbol
end

struct ProbabilityBox <: ImpreciseUQInput
    lb::AbstractVector{Real}
    ub::AbstractVector{Real}
    dist::Function
    name::Symbol
end

l = ProbabilityBox([1.75, 1.83], [1.77, 1.85], x -> Uniform(x...), :l) # length
b = Interval(0.10, 0.14, :b) # width

# h = ProbabilityBox([22.5, 0.14], [23.5, 0.16], x -> LogNormal(x...), :h) # height
h = RandomVariable(Normal(0.24, 0.01), :h) # height

μ = log(10e9^2 / sqrt(1.6e9^2 + 10e9^2))
σ = sqrt(log(1.6e9^2 / 10e9^2 + 1))
E = RandomVariable(LogNormal(μ, σ), :E) # young's modulus

μ = log(5000^2 / sqrt(400^2 + 5000^2))
σ = sqrt(log(400^2 / 5000^2 + 1))
P = RandomVariable(LogNormal(μ, σ), :P) # tip load

μ = log(600^2 / sqrt(140^2 + 600^2))
σ = sqrt(log(140^2 / 600^2 + 1))
ρ = RandomVariable(LogNormal(μ, σ), :ρ) # density

c = GaussianCopula([1 0.8; 0.8 1])
jd = JointDistribution([E, ρ], c)

inputs = [l, b, h, P, jd]

inertia = Model(df -> df.b .* df.h .^ 3 / 12, :I)

displacement = Model(
    df ->
        (df.ρ .* 9.81 .* df.b .* df.h .* df.l .^ 4) ./ (8 .* df.E .* df.I) .+
        (df.P .* df.l .^ 3) ./ (3 .* df.E .* df.I),
    :w,
)

max_displacement = 0.01

imprecise_inputs = filter(x -> isa(x, ImpreciseUQInput), inputs)
precise_inputs = filter(x -> isa(x, PreciseUQInput), inputs)

function bounds(inputs::AbstractVector{UQInput})
    lb = Float64[]
    ub = Float64[]

    for i in inputs
        if isa(i, Interval)
            push!(lb, i.lb)
            push!(ub, i.ub)
        elseif isa(i, ProbabilityBox)
            append!(lb, i.lb)
            append!(ub, i.ub)
        end
    end

    return lb, ub
end

lb, ub = bounds(inputs)

function map_to_precise_inputs(x::AbstractVector, inputs::AbstractVector{UQInput})
    precise_inputs = PreciseUQInput[]

    params = copy(x)

    for i in inputs
        if isa(i, Interval)
            par = Parameter(popfirst!(params), i.name)
            push!(precise_inputs, par)
        elseif isa(i, ProbabilityBox)
            d = length(i.lb)
            p = [popfirst!(params) for _ in 1:d]
            rv = RandomVariable(i.dist(p), i.name)
            push!(precise_inputs, rv)
        end
    end

    return precise_inputs
end

function f(x)
    inputs = [precise_inputs..., map_to_precise_inputs(x, imprecise_inputs)...]

    mc = MonteCarlo(10^6)

    mc_pf, _, _ = probability_of_failure(
        [inertia, displacement], df -> max_displacement .- df.w, inputs, mc
    )

    return mc_pf
end

x0 = (lb + ub) ./ 2

x, info = prima(f, x0; xl=lb, xu=ub)

pf_min = info.fx

x, info = prima(x -> -f(x), x0; xl=lb, xu=ub)

pf_max = -info.fx
