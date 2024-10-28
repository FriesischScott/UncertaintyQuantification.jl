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

firstorder_analytical = [0.219, 0.687, 0.00]
totaleffect_analytical = [0.3136, 0.687, 0.0946]

#==
Now, we'll create an RQMC-struct of our choice.
==#

rqmc = SobolSampling(1000, :matousek)
mc = MonteCarlo(1000)
qmc = SobolSampling(1000, :none)
#==
Then, calcualte the sobol indices:
==#

si = sobolindices(ishigami, inputs, :y, rqmc)
simc = sobolindices(ishigami, inputs, :y, mc)
siq = sobolindices(ishigami, inputs, :y, qmc)
print(reduce(+, (si.FirstOrder .- firstorder_analytical) .^ 2))
print(reduce(+, (simc.FirstOrder .- firstorder_analytical) .^ 2))
print(reduce(+, (siq.FirstOrder .- firstorder_analytical) .^ 2))

print(reduce(+, (si.TotalEffect .- totaleffect_analytical) .^ 2))
print(reduce(+, (simc.TotalEffect .- totaleffect_analytical) .^ 2))
print(reduce(+, (siq.TotalEffect .- totaleffect_analytical) .^ 2))