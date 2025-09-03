@testset "GaussianProcessDataStandardizer" begin
    single_input = RandomVariable(Normal(-1, 0.5), :x1)
    multi_input = RandomVariable.([Uniform(-2, 0), Normal(-1, 0.5), Uniform(0, 1)], [:x1, :x2, :x3])

    N = 10
    output = :y

    df_single_in = sample(single_input, N)
    df_single_in[!, output] = rand(N)

    df_multi_in = sample(multi_input, N)
    df_multi_in[!, output] = rand(N)

    single_input_names = propertynames(df_single_in[:, Not(output)])
    multi_input_names = propertynames(df_multi_in[:, Not(output)])

    @testset "IdentityTransform" begin
        dts_single = DataStandardizer(
            df_single_in, single_input_names, output, 
            InputTransform(IdentityTransform()), 
            OutputTransform(IdentityTransform())
        )

        dts_multi = DataStandardizer(
            df_multi_in, multi_input_names, output, 
            InputTransform(IdentityTransform()), 
            OutputTransform(IdentityTransform())
        )

        single_in_transformed = dts_single.fᵢ(df_single_in)
        multi_in_transformed = dts_single.fᵢ(df_single_in)
    end

    @testset "ZScoreTransform" begin
        
    end

    @testset "UnitRangeTransform" begin
        
    end

    @testset "StandardNormalTransform" begin
        #: TODO
    end
end
