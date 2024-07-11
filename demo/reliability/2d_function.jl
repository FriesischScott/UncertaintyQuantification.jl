using UncertaintyQuantification

X1 = ProbabilityBox{Normal}([Interval(-1, 2, :μ), Parameter(1, :σ)], :X1)
X2 = ProbabilityBox{Normal}([Interval(-2, 1, :μ), Parameter(2, :σ)], :X2)

model = Model(df -> 7 .- df.X1 .^ 2 .+ df.X2, :f)
inputs = [X1, X2]
simulation = MonteCarlo(10^6)
performance = df -> df.f
pf1 = probability_of_failure(model, performance, inputs, DoubleLoop(simulation))

pf2 = probability_of_failure(model, performance, inputs, IntervalMonteCarlo(5000))
