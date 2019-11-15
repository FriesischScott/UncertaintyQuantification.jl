using UncertaintyQuantification, DataFrames

l = Parameter(1.8, "l")
b = Parameter(0.12, "b")

h = RandomVariable(Normal(0.24, 0.01), "h")
E = RandomVariable(
    LogNormal(
        log(10e9^2 / sqrt(1.6e9^2 + 10e9^2)),
        sqrt(log(1.6e9^2 / 10e9^2 + 1)),
    ),
    "E",
)

P = RandomVariable(
    LogNormal(
        log(5000^2 / sqrt(400^2 + 5000^2)),
        sqrt(log(400^2 / 5000^2 + 1)),
    ),
    "P",
)
ρ = RandomVariable(
    LogNormal(log(600^2 / sqrt(140^2 + 600^2)), sqrt(log(140^2 / 600^2 + 1))),
    "ρ",
)

rvset = RandomVariableSet([ρ, E], [1 0.8; 0.8 1])

inertia = Model(df -> df.b .* df.h .^ 3 / 12, "I")

displacement = Model(
    df -> (df.ρ .* 9.81 .* df.b .* df.h .* df.l .^ 4) ./ (8 .* df.E .* df.I) .+
          (df.P .* df.l .^ 3) ./ (3 .* df.E .* df.I),
    "w",
)

inputs = [l, b, P, h, rvset]

x = mean(inputs)

g = gradient_in_standard_normal_space(
    [inertia, displacement, Model(x -> 0.01 .- x.w, "g")],
    inputs,
    x,
    :g,
)

ls = LineSampling(50, float(collect(0.5:0.5:3.5)), g)

max_displacement = 0.01

pf, samples = probability_of_failure(
    [inertia, displacement],
    df -> max_displacement .- df.w,
    inputs,
    ls,
)

println("Probability of failure: ", pf)
