@testset "GaussianProcessDataManipulation" begin
    single_input = RandomVariable(Normal(-1, 0.5), :x1)
    single_input_vector = [single_input]
    multi_input = RandomVariable.([Uniform(-2, 0), Normal(-1, 0.5), Uniform(0, 1)], [:x1, :x2, :x3])

    df_single = sample(single_input, 10)
    df_single_vector = sample(single_input_vector, 10)
    df_multi = sample(multi_input, 10)

    @testset "InputTransformer" begin
        # Check 1D input
        single_input_transformer_no = InputTransformer(df_single, names(single_input), false)
        single_input_transformer_zsc = InputTransformer(df_single, names(single_input), true)

        @test all(single_input_transformer_no(df_single) .== df_single[:, 1])
        @test all(
            single_input_transformer_zsc(df_single) .≈
            (df_single[:, 1] .- mean(df_single[:, 1])) / std(df_single[:, 1])
            )

        # Check 1D input passed in a Vector
        single_input_vector_transformer_no = InputTransformer(
            df_single_vector, names(single_input_vector), false
            )
        single_input_vector_transformer_zsc = InputTransformer(
            df_single_vector, names(single_input_vector), true
            )
        
        @test all(single_input_vector_transformer_no(df_single_vector) .== df_single_vector[:, 1])
        @test all(
            single_input_vector_transformer_zsc(df_single_vector) .≈ 
            (df_single_vector[:, 1] .- mean(df_single_vector[:, 1])) / std(df_single_vector[:, 1])
            )

        # Check ND input
        multi_input_transformer_no = InputTransformer(df_multi, names(multi_input), false)
        multi_input_transformer_zsc = InputTransformer(df_multi, names(multi_input), true)

        df_as_matrix = Matrix(df_multi)
        mean_ = mean(df_as_matrix; dims=1)
        std_ = std(df_as_matrix; dims=1)
        for (i, col) in enumerate(eachcol(df_as_matrix))
            df_as_matrix[:, i] .= (col .- mean_[1, i]) / std_[1, i]
        end

        @test all(multi_input_transformer_no(df_multi) .== Matrix(df_multi))
        @test all(multi_input_transformer_zsc(df_multi) .≈ df_as_matrix)
    end

    @testset "UQInputTransformer" begin
        # Check 1D input
        single_input_transformer_no = UQInputTransformer(single_input, false)
        single_input_transformer_sns = UQInputTransformer(single_input, true)

        df_copy_sns = copy(df_single)
        to_standard_normal_space!(single_input, df_copy_sns)

        @test all(single_input_transformer_no(df_single) .== df_single[:, 1])
        @test all(single_input_transformer_sns(df_single) .== df_copy_sns[:, 1])

        # Check 1D input passed in a Vector
        single_input_vector_transformer_no = UQInputTransformer(single_input_vector, false)
        single_input_vector_transformer_sns = UQInputTransformer(single_input_vector, true)

        df_copy_sns = copy(df_single_vector)
        to_standard_normal_space!(single_input_vector, df_copy_sns)

        @test all(single_input_transformer_no(df_single_vector) .== df_single_vector[:, 1])
        @test all(single_input_transformer_sns(df_single_vector) .== df_copy_sns[:, 1])

        # Check ND input
        multi_input_transformer_no = UQInputTransformer(multi_input, false)
        multi_input_transformer_sns = UQInputTransformer(multi_input, true)

        df_copy_sns = copy(df_multi)
        to_standard_normal_space!(multi_input, df_copy_sns)

        @test all(multi_input_transformer_no(df_multi) .== Matrix(df_multi))
        @test all(multi_input_transformer_sns(df_multi) .== Matrix(df_copy_sns))
    end

    @testset "OutputTransformer" begin
        #: TODO
    end
end
