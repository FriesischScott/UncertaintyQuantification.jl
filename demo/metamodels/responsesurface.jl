using UncertaintyQuantification
using DataFrames





x = RandomVariable.(Uniform(-π ,π), [:x1, :x2, :x3])

polynomial = Model(
    df ->
        2 .* df.x1.^2 .+ df.x1 .* df.x2 + df.x1 .* df.x3 .+ 3 .* df.x2.^2 .+ df.x2 .* df.x3 .+ 4 .* df.x3.^2 + df.x1 .+ df.x2 .+ df.x3 .+ 3,
    :y,
)

training_data = sample(x, 100)

evaluate!(polynomial, training_data)

rs = ResponseSurface(training_data, :y, 2)

test_data = sample(x, 1000)
evaluate!(rs, test_data)

p_data = test_data[:, Not(:y)]
evaluate!(polynomial, p_data)

mse = mean((p_data.y .- test_data.y) .^ 2)
println("MSE is:  $mse")



