@testset "Parameter" begin

    π = Parameter(3.14, "π")

    @test typeof(π) == Parameter
    @test π.value == 3.14
    @test π.name == "π"

end
