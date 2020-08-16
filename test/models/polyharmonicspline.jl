@testset "PolyharmonicSpline" begin
    A = [1 2; 1 3; 3 6; 4 3; 5 7; 5 2; 6 8; 6 9; 3 7]
    f = [2; 3; 18; 12; 35; 10; 48; 54; 21]
    k = 1
    order = df -> [df.a df.b]

    input = DataFrame(a = 4, b = 3)

    polyspline = PolyharmonicSpline(A, f, k, :fkt, order)

    B = Array{Float64}(undef, 1)
    B[1] = 12

    @test isa(polyspline, PolyharmonicSpline)

    B = Array{Float64}(undef, 1)
    B[1] = 12

    evaluate!(polyspline, input)

    A = input.fkt
    A = round.(A, digits=3)

    @test A == B
end
