using UncertaintyQuantification

form = FORM()

x1 = RandomVariable(Normal(200, 20), :x1)
x2 = RandomVariable(Normal(150, 10), :x2)

inputs = [x1, x2]

@time pf, β, dp = probability_of_failure(df -> df.x1 .- df.x2, inputs, form)
