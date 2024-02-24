@testset "Interval propagation" begin

    X1 = ProbabilityBox{Normal}([Interval(-1, 2, :μ), Parameter(1, :σ)], :X1) 
    X2 = ProbabilityBox{Normal}([Interval(-2, 1, :μ), Parameter(2, :σ)], :X2) 
    X3 = RandomVariable(Normal(0, 1), :X3)
    X4 = Parameter(5, :X4)

    inputs = [X1, X2, X3, X4]
    models = Model(df ->  7 .− df.X1.^2 .+ df.X2 .+ df.X3 .+ df.X4, :g)
    performance = df -> df.g

    df = sample(inputs, 1000)
    evaluate!(models, df)

    @test eltype(df.g) == Interval


end