using UncertaintyQuantification, FiniteDifferences

l = Parameter(1.8, "l")
b = Parameter(0.12, "b")

h = RandomVariable(Normal(0.24, 0.01), "h")
E = RandomVariable(
    LogNormal(
        log(10e9^2 / sqrt(1.6e9 + 10e9^2)),
        sqrt(log(1.6e9 / 10e9^2 + 1)),
    ),
    "E",
)

P = RandomVariable(
    LogNormal(log(5000^2 / sqrt(400 + 5000^2)), sqrt(log(400 / 5000^2 + 1))),
    "P",
)
ρ = RandomVariable(
    LogNormal(log(600^2 / sqrt(140 + 600^2)), sqrt(log(140 / 600^2 + 1))),
    "ρ",
)

rvset = RandomVariableSet([P ρ], [1 0.8; 0.8 1])

inertia = Model(df -> df.b .* df.h .^ 3 / 12, "I")

displacement = Model(
    df -> (df.ρ .* 9.81 .* df.b .* df.h .* df.l .^ 4) ./
          (8 .* df.E .* df.I) .+ (df.P .* df.l .^ 3) ./ (3 .* df.E .* df.I),
    "w",
)

x = map(mean, [h E P ρ])

g = gradient([inertia, displacement], [l, b, h, E, rvset], x, :w)

println("Gradient: $g")
