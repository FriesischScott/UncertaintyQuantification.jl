@testset "GaussianProcessRegression" begin

    @testset "OneDimensionalInput" begin
        # ---------------------------------------------------
        # Test construction from DataFrame
        x = collect(range(0, stop=5, length=10))
        y = sin.(x)
        data = DataFrame(:x => x, :y => y)

        σ² = 1e-5
        kernel = SqExponentialKernel()

        gp = GP(0.0, kernel)
        gpr = GaussianProcess(
            gp,
            data,
            :y;
            optimization=NoHyperparameterOptimization()
        )

        test_data = DataFrame(:x => x)
        evaluate!(gpr, test_data)
        mean_and_var!(gpr, test_data)
        
        # evaluate! returns mean as standard
        @test all(test_data[!, :y] .== test_data[!, :y_mean])
        # outputs at trainingset should be very close
        @test all(isapprox.(test_data[!, :y], y; atol=1e-14))
        # variance should be very close to zero as we did not use observation noise
        @test all(isapprox.(test_data[!, :y_var], 0.0; atol=1e-14))

        noisy_gp = with_gaussian_noise(GP(0.0, kernel), σ²)
        noisy_gpr = GaussianProcess(
            noisy_gp,
            data,
            :y;
            optimization=NoHyperparameterOptimization()
        )

        test_data = DataFrame(:x => x)
        evaluate!(noisy_gpr, test_data)
        mean_and_var!(noisy_gpr, test_data)

        # evaluate! returns mean as standard
        @test all(test_data[!, :y] .== test_data[!, :y_mean])
        # check if prediction variance is within 5% deviation from prescribed noise
        @test all(abs.(test_data[!, :y_var] .- σ²) .< 0.05σ²)

        # ---------------------------------------------------
        # Test construction from UQInput + UQModel
        x = RandomVariable(Uniform(0, 5), :x)
        model = Model(
            df -> sin.(df.x), :y
        )
        design = LatinHypercubeSampling(10)

        σ² = 1e-5
        kernel = SqExponentialKernel()

        gp = GP(0.0, kernel)
        gpr = GaussianProcess(
            gp,
            x,
            model,
            :y,
            design;
            optimization=NoHyperparameterOptimization()
        )

        test_data = sample(x, design)
        evaluate!(gpr, test_data)
        mean_and_var!(gpr, test_data)
        
        # evaluate! returns mean as standard
        @test all(test_data[!, :y] .== test_data[!, :y_mean])

        noisy_gp = with_gaussian_noise(GP(0.0, kernel), σ²)
        noisy_gpr = GaussianProcess(
            noisy_gp,
            x,
            model,
            :y,
            design;
            optimization=NoHyperparameterOptimization()
        )

        test_data = sample(x, design)
        evaluate!(gpr, test_data)
        mean_and_var!(gpr, test_data)

        # evaluate! returns mean as standard
        @test all(test_data[!, :y] .== test_data[!, :y_mean])
    end
    @testset "TwoDimensionalInput" begin
        # ---------------------------------------------------
        # Test construction from DataFrame
        x = [collect(range(0, stop=5, length=10)), collect(range(0, stop=5, length=10))]
        y = sin.(x[1]) + cos.(x[2])
        data = DataFrame(:x1 => x[1], :x2 => x[2], :y => y)

        σ² = 1e-5
        kernel = SqExponentialKernel()

        gp = GP(0.0, kernel)
        gpr = GaussianProcess(
            gp,
            data,
            :y;
            optimization=NoHyperparameterOptimization()
        )

        test_data = DataFrame(:x1 => x[1], :x2 => x[2])
        evaluate!(gpr, test_data)
        mean_and_var!(gpr, test_data)
        
        # evaluate! returns mean as standard
        @test all(test_data[!, :y] .== test_data[!, :y_mean])
        # outputs at trainingset should be very close
        @test all(isapprox.(test_data[!, :y], y; atol=1e-14))
        # variance should be very close to zero as we did not use observation noise
        @test all(isapprox.(test_data[!, :y_var], 0.0; atol=1e-14))

        noisy_gp = with_gaussian_noise(GP(0.0, kernel), σ²)
        noisy_gpr = GaussianProcess(
            noisy_gp,
            data,
            :y;
            optimization=NoHyperparameterOptimization()
        )

        test_data = DataFrame(:x1 => x[1], :x2 => x[2])
        evaluate!(noisy_gpr, test_data)
        mean_and_var!(noisy_gpr, test_data)

        # evaluate! returns mean as standard
        @test all(test_data[!, :y] .== test_data[!, :y_mean])
        # check if prediction variance is within 5% deviation from prescribed noise
        @test all(abs.(test_data[!, :y_var] .- σ²) .< 0.05σ²)

        # ---------------------------------------------------
        # Test construction from UQInput + UQModel
        x = RandomVariable.([Uniform(0, 5), Uniform(0, 5)], [:x1, :x2])
        model = Model(
            df -> sin.(df.x1) + cos.(df.x2), :y
        )
        design = LatinHypercubeSampling(10)

        σ² = 1e-5
        kernel = SqExponentialKernel()

        gp = GP(0.0, kernel)
        gpr = GaussianProcess(
            gp,
            x,
            model,
            :y,
            design;
            optimization=NoHyperparameterOptimization()
        )

        test_data = sample(x, design)
        evaluate!(gpr, test_data)
        mean_and_var!(gpr, test_data)
        
        # evaluate! returns mean as standard
        @test all(test_data[!, :y] .== test_data[!, :y_mean])

        noisy_gp = with_gaussian_noise(GP(0.0, kernel), σ²)
        noisy_gpr = GaussianProcess(
            noisy_gp,
            x,
            model,
            :y,
            design;
            optimization=NoHyperparameterOptimization()
        )

        test_data = sample(x, design)
        evaluate!(gpr, test_data)
        mean_and_var!(gpr, test_data)

        # evaluate! returns mean as standard
        @test all(test_data[!, :y] .== test_data[!, :y_mean])
    end
end