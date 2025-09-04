function make_standardizer(
    dataframe, 
    input_names, 
    output,
    it,
    ot
)
    UncertaintyQuantification.DataStandardizer(
        dataframe, input_names, output, 
        UncertaintyQuantification.InputTransform(it),
        UncertaintyQuantification.OutputTransform(ot),
    )
end

@testset "GaussianProcessDataStandardizer" begin
    transforms = [
        UncertaintyQuantification.IdentityTransform(),
        UncertaintyQuantification.ZScoreTransform(),
        UncertaintyQuantification.UnitRangeTransform(),
        UncertaintyQuantification.StandardNormalTransform()
    ]

    N = 10
    output = :y

    single_input = RandomVariable(Normal(-1, 0.5), :x1)
    df_single = sample(single_input, N)
    df_single[!, output] = rand(N)

    multi_input = RandomVariable.([Uniform(-2, 0), Normal(-1, 0.5), Uniform(0, 1)], [:x1, :x2, :x3])
    df_multi = sample(multi_input, N)
    df_multi[!, output] = df_single[!, output]
    
    names_single = propertynames(df_single[:, Not(output)])
    names_multi = propertynames(df_multi[:, Not(output)])

    for transform in transforms
        @testset "$(nameof(typeof(transform)))" begin
            for (testname, df, inputs, names) in [
                    ("single input", df_single, single_input, names_single),
                    ("multi input",  df_multi,  multi_input, names_multi)
                ]
                @testset "$testname" begin
                    if isa(transform, UncertaintyQuantification.StandardNormalTransform)
                        @test_throws ArgumentError UncertaintyQuantification.DataStandardizer(
                            df, inputs, output,
                            UncertaintyQuantification.InputTransform(transform),
                            UncertaintyQuantification.OutputTransform(transform)
                        )
                        # Test output shapes for identity output transform to check StandardNormalTransform for inputs
                        dts = UncertaintyQuantification.DataStandardizer(
                            df, inputs, output,
                            UncertaintyQuantification.InputTransform(transform),
                            UncertaintyQuantification.OutputTransform(
                                UncertaintyQuantification.IdentityTransform()
                                )
                        )
                        Xin = dts.fᵢ(df)
                        if testname == "single input"
                            # input gets transformed to a Vector
                            @test isa(Xin, Vector)
                        else
                            # input gets transformed to RowVecs
                            @test isa(Xin, RowVecs)
                        end
                    else
                        dts = UncertaintyQuantification.DataStandardizer(
                            df, inputs, output,
                            UncertaintyQuantification.InputTransform(transform),
                            UncertaintyQuantification.OutputTransform(transform),
                        )

                        Xin = dts.fᵢ(df)
                        Xout = dts.fₒ(df)
                        Yout = dts.fₒ⁻¹(Xout)
                        var_Yout = dts.var_fₒ⁻¹(Xout)
                        if testname == "single input"
                            # input gets transformed to a Vector
                            @test isa(Xin, Vector)
                            if isa(transform, UncertaintyQuantification.IdentityTransform)
                                # Test input scaling
                                @test all(Xin .== df[!, only(names)])
                                # Test output scaling
                                @test all(Yout .== Xout)
                                # Test output inverse transform
                                @test all(df[!, output] .== Yout)
                                # Test output inverse transform for variance
                                @test all(df[!, output] .== var_Yout)

                            elseif isa(transform, UncertaintyQuantification.ZScoreTransform)
                                # Test input scaling
                                μ = mean(df[!, only(names)])
                                σ = std(df[!, only(names)]) 
                                Min = (df[!, only(names)] .- μ) ./ σ
                                @test all(Xin .≈ Min)

                                # Test output scaling
                                μ = mean(df[!, output])
                                σ = std(df[!, output]) 
                                Mout = (df[!, output] .- μ) ./ σ                                    
                                @test all(Mout .≈ Xout)
                                # Test output inverse transform
                                @test all(df[!, output] .≈ Yout)
                                # Test output inverse transform for variance
                                @test all(σ^2 * Xout .≈ var_Yout)

                            elseif isa(transform, UncertaintyQuantification.UnitRangeTransform)
                                # Test input scaling
                                tmin, tmax = extrema(df[!, only(names)])
                                shift = tmin
                                scale = 1 / (tmax - tmin)
                                Min = (df[!, only(names)] .- shift) * scale
                                @test all(Xin .≈ Min)

                                # Test output scaling
                                tmin, tmax = extrema(df[!, output])
                                shift = tmin
                                scale = 1 / (tmax - tmin)
                                Mout = (df[!, output] .- shift) * scale                                   
                                @test all(Mout .≈ Xout)
                                # Test output inverse transform
                                @test all(df[!, output] .≈ Yout)
                                # Test output inverse transform for variance
                                @test all(scale^2 * Xout .≈ var_Yout)

                            end
                        else
                            # input gets transformed to RowVecs
                            @test isa(Xin, RowVecs)
                            if isa(transform, UncertaintyQuantification.IdentityTransform)
                                # Test input scaling
                                Min = mapreduce(rv -> rv', vcat, Xin)
                                @test all(Min .== Matrix(df[!, names]))
                                # Test output scaling
                                @test all(Yout .== Xout)
                                # Test output inverse transform
                                @test all(df[!, output] .== Yout)
                                # Test output inverse transform for variance
                                @test all(df[!, output] .== var_Yout)

                            elseif isa(transform, UncertaintyQuantification.ZScoreTransform)
                                # Test input scaling
                                Xin = mapreduce(rv -> rv', vcat, Xin)
                                μ = mean(Matrix(df[!, names]), dims=1)
                                σ = std(Matrix(df[!, names]), dims=1) 
                                Min = (Matrix(df[!, names]) .- μ) ./ σ
                                @test all(Xin .≈ Min)

                                # Test output scaling
                                μ = mean(df[!, output])
                                σ = std(df[!, output]) 
                                Mout = (df[!, output] .- μ) ./ σ                                    
                                @test all(Mout .≈ Xout)
                                # Test output inverse transform
                                @test all(df[!, output] .≈ Yout)
                                # Test output inverse transform for variance
                                @test all(σ^2 * Xout .≈ var_Yout)

                            elseif isa(transform, UncertaintyQuantification.UnitRangeTransform)
                                # Test input scaling
                                Xin = mapreduce(rv -> rv', vcat, Xin)
                                extrs = extrema(Matrix(df[!, names]), dims=1)
                                shift = map(t -> t[1], extrs[1, :])
                                scale = map(t -> 1 / (t[2] - t[1]), extrs[1, :])
                                Min = (Matrix(df[!, names]) .- shift') .* scale'
                                @test all(Xin .≈ Min)

                                # Test output scaling
                                tmin, tmax = extrema(df[!, output])
                                shift = tmin
                                scale = 1 / (tmax - tmin)
                                Mout = (df[!, output] .- shift) * scale                                   
                                @test all(Mout .≈ Xout)
                                # Test output inverse transform
                                @test all(df[!, output] .≈ Yout)
                                # Test output inverse transform for variance
                                @test all(scale^2 * Xout .≈ var_Yout)

                            end
                        end
                    end
                end
            end
        end
    end
end
