using UncertaintyQuantification

inputs = Vector{RandomVariable}(undef, 20)

X = [1,1.1,0.9,1.2,0.8,0.7,1.3,1.4,0.6,0.95]

for i in 1:10
    inputs[i] = RandomVariable(Normal(0, X[i]), Symbol("A$i"))
    inputs[i+10] = RandomVariable(Normal(0, X[i]), Symbol("B$i"))
end


bfkt = Model(df -> df.A1 .* df.A6 + df.A2 .* df.A7 + df.A3 .* df.A8
+ df.A4 .* df.A9 + df.A5 .* df.A10, :fkt)


mc = MonteCarlo(10000)

si = soindex([bfkt], inputs, mc)
