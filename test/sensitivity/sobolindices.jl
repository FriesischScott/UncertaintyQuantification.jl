@testset "sobolindices" begin

inputs = Vector{RandomVariable}(undef, 10)

X = [1,1.1,0.9,1.2,0.8,0.7,1.3,1.4,0.6,0.95]

for i in 1:5
    inputs[i] = RandomVariable(Normal(0, X[i]), Symbol("x$i"))
    inputs[i+5] = RandomVariable(Normal(0, X[i+5]), Symbol("w$i"))
end


bfkt = Model(df -> df.x1 .* df.w1 + df.x2 .* df.w2 + df.x3 .* df.w3
+ df.x4 .* df.w4 + df.x5 .* df.w5, :fkt)


mc = MonteCarlo(10000)

si = sobolindices([bfkt], inputs, mc)

round_si = round.(si, digits = 1)

round_analytical = [0.1; 0.4; 0.3; 0.1; 0.1; 0.1; 0.4; 0.3; 0.1; 0.1]

@test round_si == round_analytical

end
