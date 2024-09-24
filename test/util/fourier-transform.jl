@testset "Fourier transform" begin

    t = 0:0.01:10
    x = sin.(2*π*t)
    
    X = dft(x)

    x_hat = idft(X)
    x_approx = real(x_hat)

    @test x≈x_approx

end
 