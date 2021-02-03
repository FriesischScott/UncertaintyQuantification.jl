using UncertaintyQuantification, Plots

x = RandomVariable(Normal(), :x)
y = RandomVariable(Normal(), :y)

model = Model(df -> 0.35 * sqrt.(df.x.^2 .+ df.y.^2), :z)

subset = UncertaintyQuantification.SubSetSimulation(1000, 0.1, 10, 0.5)

pf, samples = probability_of_failure(
    [model],
    df -> 1 .- df.z,
    [x, y],
    subset)

println("Probability of failure: $pf")

p = plot()

for level ∈ unique(samples.level)
    data = filter(:level => l -> l == level, samples)
    to_standard_normal_space!([x, y], data)
    scatter!(p, data.x, data.y, label="l_$level")
end

display(p)