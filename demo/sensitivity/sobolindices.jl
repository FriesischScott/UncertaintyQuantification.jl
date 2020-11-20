using UncertaintyQuantification

σx = [1, 1.1, 0.9, 1.2, 0.8]
σω = [0.7, 1.3, 1.4, 0.6, 0.95]

X = RandomVariable.(Normal.(0, σx), [:X1, :X2, :X3, :X4, :X5])
ω = RandomVariable.(Normal.(0, σω), [:ω1, :ω2, :ω3, :ω4, :ω5])

B = Model(df -> df.X1 .* df.ω1
    + df.X2 .* df.ω2
    + df.X3 .* df.ω3
    + df.X4 .* df.ω4
    + df.X5 .* df.ω5, :B)

mc = MonteCarlo(10000)

s_mc, st_mc = sobolindices([B], [X; ω], :B, mc)

# Compare with analytical solution
VB = sum(σx.^2 .* σω.^2)
st_analytical = (σx.^2 .* σω.^2) / VB
# Monte Carlo
@assert isapprox(s_mc.X1, 0, atol=0.1)
@assert isapprox(s_mc.X2, 0, atol=0.1)
@assert isapprox(s_mc.X3, 0, atol=0.1)
@assert isapprox(s_mc.X4, 0, atol=0.1)
@assert isapprox(s_mc.X5, 0, atol=0.1)

@assert isapprox(st_mc.X1, st_analytical[1], rtol=0.1)
@assert isapprox(st_mc.X2, st_analytical[2], rtol=0.1)
@assert isapprox(st_mc.X3, st_analytical[3], rtol=0.1)
@assert isapprox(st_mc.X4, st_analytical[4], rtol=0.1)
@assert isapprox(st_mc.X5, st_analytical[5], rtol=0.1)