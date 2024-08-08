#==
## Quasi Monte Carlo 

## Example

In this example, we will perform a semsitivity analysis using Quasi Monte Carlo.
Let's consider the ishigami function which is defined as $f(x_{1}, x_{2}, x_{3}) = sin(x_{1}) +  a * sin(x_{2})^2 + b * x_{3}^4 * sin(x_{1}) $
==#
using UncertaintyQuantification

x = RandomVariable.(Uniform(-π, π), [:x1, :x2, :x3])
a = Parameter(7, :a)
b = Parameter(0.05, :b)

inputs = [x; a; b]

ishigami = Model(
    df -> sin.(df.x1) .+ df.a .* sin.(df.x2) .^ 2 .+ df.b .* df.x3 .^ 4 .* sin.(df.x1), :y
)

#==
Now, we'll create an RQMC-struct of our choice.
==#

rqmc = SobolSampling(1000, :owenscramble)

#==
Then, calcualte the sobol indices:
==#

si = sobolindices(ishigami, inputs, :y, rqmc)
