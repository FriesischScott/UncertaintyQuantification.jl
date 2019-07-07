include("../src/OpenCossan.jl")

using Main.OpenCossan
using Distributions, DataFrames

import LinearAlgebra.I

nmc = 10000

l = 1.8
b = 0.12
max_displacement = 0.01

h = Normal(0.24, 0.01)
P = LogNormal(log(5000^2/sqrt(400 + 5000^2)), sqrt(log(400/5000^2 + 1)))
rho = LogNormal(log(600^2/sqrt(140 + 600^2)), sqrt(log(140/600^2 + 1)))
E = LogNormal(log(10e9^2/sqrt(1.6e9 + 10e9^2)), sqrt(log(1.6e9/10e9^2 + 1)))

corr = Matrix{Float64}(I, 4, 4)
corr[3,2] = 0.8
corr[2,3] = 0.8

rvset = RandomVariableSet([h P rho E], ["h" "P" "rho" "E"], corr)

samples = DataFrame(l = fill(l, nmc), b = fill(b, nmc))

samples = [samples rand(rvset, nmc)]

function inertia(df::DataFrame)
    df.b .* df.h .^ 3 / 12
end

samples.I = inertia(samples)

function displacement(df::DataFrame)
    (df.rho .* 9.81 .* df.b .* df.h .* df.l .^ 4) ./ (8 .* df.E .* df.I) .+ (df.P .* df.l .^ 3) ./ (3 .* df.E .* df.I)
end

samples.w = displacement(samples)

performance = -samples.w .+ max_displacement

pf = sum(performance .< 0) / nmc

println("Probability of failure: ", pf)


