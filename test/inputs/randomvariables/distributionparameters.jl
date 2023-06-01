@testset "Distribution Parameters" begin
    @test [distribution_parameters(80, 15, LogNormal)...] ≈
        [4.364750443920552, 0.18588270900398385]
end
