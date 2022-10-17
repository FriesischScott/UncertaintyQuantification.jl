using UncertaintyQuantification
using DataFrames

x = RandomVariable.(Uniform(-10 , 10), [:x1, :x2])

inputs = [x[1]; x[2]]

polynomial = Model(
    df ->
    2 *(df.x1).^2 .+ 0.5 * (df.x2).^2 .+ df.x1 .* df.x2 .+ df.x1 .+ 8 * df.x2 .+ 5,
    :y,
)

training_data = sample(inputs, 100)

evaluate!(polynomial, training_data)

rs = ResponseSurface(training_data, :y, 2, 3)

test_data = sample(x, 1000)
evaluate!(rs, test_data)

p_data = test_data[:, Not(:y)]
evaluate!(polynomial, p_data)

mse = mean((p_data.y .- test_data.y) .^ 2)
println("MSE is:  $mse")



