@testset "sobolindices" begin

    σx = [1, 1.1, 0.9, 1.2, 0.8]
    σω = [0.7, 1.3, 1.4, 0.6, 0.95]

    X = RandomVariable.(Normal.(0, σx), [:X1, :X2, :X3, :X4, :X5])
    ω = RandomVariable.(Normal.(0, σω), [:ω1, :ω2, :ω3, :ω4, :ω5])

    B = Model(df -> df.X1 .* df.ω1
    + df.X2 .* df.ω2
    + df.X3 .* df.ω3
    + df.X4 .* df.ω4
    + df.X5 .* df.ω5, :B)

    si = sobolindices([B], [X; ω], :B, MonteCarlo(5000))

    @test round.([si...], digits=1) == [0.1; 0.4; 0.3; 0.1; 0.1; 0.1; 0.4; 0.3; 0.1; 0.1]

end