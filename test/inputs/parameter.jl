@testset "Parameter" begin

    π = Parameter(3.14, :π)

    @test isa(π, Parameter)
    @test π.value == 3.14
    @test π.name == :π

    @test sample(π, 5) == DataFrame(π = [3.14, 3.14, 3.14, 3.14, 3.14])
end
