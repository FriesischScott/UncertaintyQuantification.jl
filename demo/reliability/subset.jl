using UncertaintyQuantification

x = RandomVariable(Normal(), :x)
y = RandomVariable(Normal(), :y)

model = Model(df -> 0.35 * sqrt.(df.x.^2 .+ df.y.^2), :z)

subset = UncertaintyQuantification.SubSetSimulation(1000, 0.1, 10, 0.5)

pf_mc, _ = probability_of_failure(
    [model],
    df -> 1 .- df.z,
    [x, y],
    MonteCarlo(1e7))

pf_subset = probability_of_failure(
    [model],
    df -> 1 .- df.z,
    [x, y],
    subset)

@show pf_mc
@show pf_subset