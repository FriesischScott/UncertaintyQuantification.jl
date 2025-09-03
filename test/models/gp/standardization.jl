@testset "GaussianProcessDataStandardizer" begin
    single_input = RandomVariable(Normal(-1, 0.5), :x1)
    multi_input = RandomVariable.([Uniform(-2, 0), Normal(-1, 0.5), Uniform(0, 1)], [:x1, :x2, :x3])

    N = 10
    output = :y

    df_single_in = sample(single_input, N)
    df_single_in[!, output] = rand(N)

    df_multi_in = sample(multi_input, N)
    df_multi_in[!, output] = df_single_in[!, output]

    single_input_names = propertynames(df_single_in[:, Not(output)])
    multi_input_names = propertynames(df_multi_in[:, Not(output)])

    @testset "IdentityTransform" begin
        # 1D input
        dts_single = UncertaintyQuantification.DataStandardizer(
            df_single_in, single_input_names, output, 
            UncertaintyQuantification.InputTransform(
                UncertaintyQuantification.IdentityTransform()
                ), 
            UncertaintyQuantification.OutputTransform(
                UncertaintyQuantification.IdentityTransform()
                )
        )
        single_in_transformed = dts_single.fᵢ(df_single_in)
        
        @test isa(single_in_transformed, Vector)
        @test all(single_in_transformed .== df_single_in[:, only(single_input_names)])

        # ND input
        dts_multi = UncertaintyQuantification.DataStandardizer(
            df_multi_in, multi_input_names, output, 
            UncertaintyQuantification.InputTransform(
                UncertaintyQuantification.IdentityTransform()
                ), 
            UncertaintyQuantification.OutputTransform(
                UncertaintyQuantification.IdentityTransform()
                )
        )
        multi_in_transformed = dts_multi.fᵢ(df_multi_in)

        @test isa(multi_in_transformed, RowVecs)

        RowVecsMatrix = mapreduce(rv -> rv', vcat, multi_in_transformed)
        @test all(RowVecsMatrix .== Matrix(df_multi_in[:, multi_input_names]))

        # Output
        out_transformed = dts_single.fₒ(df_single_in)

        # single input and multi input related output transforms should do the same thing
        @test all(out_transformed .== dts_multi.fₒ(df_multi_in))

        # Mean and latent function samples just get identity transformed
        @test all(df_single_in[!, output] .== dts_single.fₒ⁻¹(out_transformed))
        @test all(df_multi_in[!, output] .== dts_multi.fₒ⁻¹(out_transformed))

        # Variance also just get identity transformed
        @test all(df_single_in[!, output] .== dts_single.var_fₒ⁻¹(out_transformed))
        @test all(df_multi_in[!, output] .== dts_multi.var_fₒ⁻¹(out_transformed))
    end

    @testset "ZScoreTransform" begin
        # 1D input
        dts_single = UncertaintyQuantification.DataStandardizer(
            df_single_in, single_input_names, output, 
            UncertaintyQuantification.InputTransform(
                UncertaintyQuantification.ZScoreTransform()
                ), 
            UncertaintyQuantification.OutputTransform(
                UncertaintyQuantification.ZScoreTransform()
                )
        )
        single_in_transformed = dts_single.fᵢ(df_single_in)
        
        @test isa(single_in_transformed, Vector)

        μ = mean(df_single_in[:, only(single_input_names)])
        σ = std(df_single_in[:, only(single_input_names)])
        manually_scaled = (df_single_in[:, only(single_input_names)] .- μ) ./ σ
        @test all(single_in_transformed .≈ manually_scaled)

        # ND input
        dts_multi = UncertaintyQuantification.DataStandardizer(
            df_multi_in, multi_input_names, output, 
            UncertaintyQuantification.InputTransform(
                UncertaintyQuantification.ZScoreTransform()
                ), 
            UncertaintyQuantification.OutputTransform(
                UncertaintyQuantification.ZScoreTransform()
                )
        )
        multi_in_transformed = dts_multi.fᵢ(df_multi_in)

        @test isa(multi_in_transformed, RowVecs)

        RowVecsMatrix = mapreduce(rv -> rv', vcat, multi_in_transformed)
        μ = mean(Matrix(df_multi_in[:, multi_input_names]), dims=1)
        σ = std(Matrix(df_multi_in[:, multi_input_names]), dims=1)

        manually_scaled = (Matrix(df_multi_in[:, multi_input_names]) .- μ) ./ σ
        @test all(RowVecsMatrix .≈ manually_scaled)

        # Output
        out_transformed = dts_single.fₒ(df_single_in)

        # single input and multi input related output transforms should do the same thing
        @test all(out_transformed .== dts_multi.fₒ(df_multi_in))

        # Mean and latent function samples get rescaled and shifted
        @test all(df_single_in[!, output] .≈ dts_single.fₒ⁻¹(out_transformed))
        @test all(df_multi_in[!, output] .≈ dts_multi.fₒ⁻¹(out_transformed))

        # Variance gets multiplied by squared standard deviation used in ZScoreTransform
        # Note: This usually gets applied to the GP posterior variance, here
        # we just scale out_transformed back to verify it does the right thing
        σ = std(df_single_in[:, output])

        @test all(out_transformed * σ^2 .≈ dts_single.var_fₒ⁻¹(out_transformed))
        @test all(out_transformed * σ^2 .≈ dts_multi.var_fₒ⁻¹(out_transformed))
    end

    @testset "UnitRangeTransform" begin
        # 1D input
        dts_single = UncertaintyQuantification.DataStandardizer(
            df_single_in, single_input_names, output, 
            UncertaintyQuantification.InputTransform(
                UncertaintyQuantification.UnitRangeTransform()
                ), 
            UncertaintyQuantification.OutputTransform(
                UncertaintyQuantification.UnitRangeTransform()
                )
        )
        single_in_transformed = dts_single.fᵢ(df_single_in)
        
        @test isa(single_in_transformed, Vector)

        tmin, tmax = extrema(df_single_in[:, only(single_input_names)])
        scale = 1 / (tmax - tmin)
        manually_scaled = (df_single_in[:, only(single_input_names)] .- tmin) * scale
        @test all(single_in_transformed .≈ manually_scaled)

        # ND input
        dts_multi = UncertaintyQuantification.DataStandardizer(
            df_multi_in, multi_input_names, output, 
            UncertaintyQuantification.InputTransform(
                UncertaintyQuantification.UnitRangeTransform()
                ), 
            UncertaintyQuantification.OutputTransform(
                UncertaintyQuantification.UnitRangeTransform()
                )
        )
        multi_in_transformed = dts_multi.fᵢ(df_multi_in)

        @test isa(multi_in_transformed, RowVecs)

        RowVecsMatrix = mapreduce(rv -> rv', vcat, multi_in_transformed)
        mins_maxs = extrema(df_single_in[:, only(single_input_names)], dims=1)
        tmin = map(t -> t[1], mins_maxs[1, :])
        scale = map(t -> 1 / (t[2] - t[1]), mins_maxs[1, :])

        manually_scaled = (Matrix(df_multi_in[:, multi_input_names]) .- tmin) .* scale
        print(RowVecsMatrix)
        print(manually_scaled)
        @test all(RowVecsMatrix .≈ manually_scaled)

        # Output
        out_transformed = dts_single.fₒ(df_single_in)

        # single input and multi input related output transforms should do the same thing
        @test all(out_transformed .== dts_multi.fₒ(df_multi_in))

        # Mean and latent function samples get rescaled and shifted
        @test all(df_single_in[!, output] .≈ dts_single.fₒ⁻¹(out_transformed))
        @test all(df_multi_in[!, output] .≈ dts_multi.fₒ⁻¹(out_transformed))

        # Variance gets multiplied by squared scale used in UnitRangeTransform
        # Note: This usually gets applied to the GP posterior variance, here
        # we just scale out_transformed back to verify it does the right thing
        tmin, tmax = extrema(df_single_in[:, output])
        scale = 1 / (tmax - tmin)

        print(out_transformed * (1/scale)^2)
        print(dts_single.var_fₒ⁻¹(out_transformed))
        @test all(out_transformed * (1/scale)^2 .≈ dts_single.var_fₒ⁻¹(out_transformed))
        @test all(out_transformed * (1/scale)^2 .≈ dts_multi.var_fₒ⁻¹(out_transformed))
    end

    @testset "StandardNormalTransform" begin
        #: TODO
    end
end
