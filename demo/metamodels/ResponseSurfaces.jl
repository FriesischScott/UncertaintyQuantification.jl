using UncertaintyQuantification


x = RandomVariable.(Uniform(0, 1), [:x1, :x2])
y = RandomVariable.(Uniform(0.3, 0.7), :y)

inputs = [x; y]
training_data = sample(inputs, 100)

rs = ResponseSurface(training_data, :y, 2, 3)

test_data = sample(x, 1000)

evaluate!(rs, test_data)