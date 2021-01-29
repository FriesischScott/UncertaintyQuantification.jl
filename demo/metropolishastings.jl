using UncertaintyQuantification

h = RandomVariable(Normal(0.24, 0.01), :h)
E = RandomVariable(
    LogNormal(
        log(10e9^2 / sqrt(1.6e9^2 + 10e9^2)),
        sqrt(log(1.6e9^2 / 10e9^2 + 1)),
    ),
    :E,
)

P = RandomVariable(
    LogNormal(log(5000^2 / sqrt(400^2 + 5000^2)), sqrt(log(400^2 / 5000^2 + 1))),
    :P,
)
ρ = RandomVariable(
    LogNormal(log(600^2 / sqrt(140^2 + 600^2)), sqrt(log(140^2 / 600^2 + 1))),
    :ρ,
)

c = GaussianCopula([1 0.8; 0.8 1])
jd = JointDistribution([E ρ], c)

inputs = [h, P, jd]

mc = assemblechains(inputs, 0.1, sample(inputs, 10))

mc_mh = metropolishastings(mc, 1)