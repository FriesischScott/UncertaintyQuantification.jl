@testset "Design of Experiments" begin
    @testset "TwoLevelFactorial" begin
        tlf = TwoLevelFactorial()

        @test doe_samples(tlf, 2) == [
            0.0 0.0
            1.0 0.0
            0.0 1.0
            1.0 1.0
        ]

        @test doe_samples(tlf, 3) == [
            0.0 0.0 0.0
            1.0 0.0 0.0
            0.0 1.0 0.0
            1.0 1.0 0.0
            0.0 0.0 1.0
            1.0 0.0 1.0
            0.0 1.0 1.0
            1.0 1.0 1.0
        ]

        @test doe_samples(tlf, 4) == [
            0.0 0.0 0.0 0.0
            1.0 0.0 0.0 0.0
            0.0 1.0 0.0 0.0
            1.0 1.0 0.0 0.0
            0.0 0.0 1.0 0.0
            1.0 0.0 1.0 0.0
            0.0 1.0 1.0 0.0
            1.0 1.0 1.0 0.0
            0.0 0.0 0.0 1.0
            1.0 0.0 0.0 1.0
            0.0 1.0 0.0 1.0
            1.0 1.0 0.0 1.0
            0.0 0.0 1.0 1.0
            1.0 0.0 1.0 1.0
            0.0 1.0 1.0 1.0
            1.0 1.0 1.0 1.0
        ]
    end

    @testset "FullFactorial" begin
        ff = FullFactorial([3, 4, 2])
        @test ff.levels == [3, 4, 2]
        @test doe_samples(ff, 3) == [
            0.0 0.0 0.0
            0.5 0.0 0.0
            1.0 0.0 0.0
            0.0 0.3333333333333333 0.0
            0.5 0.3333333333333333 0.0
            1.0 0.3333333333333333 0.0
            0.0 0.6666666666666666 0.0
            0.5 0.6666666666666666 0.0
            1.0 0.6666666666666666 0.0
            0.0 1.0 0.0
            0.5 1.0 0.0
            1.0 1.0 0.0
            0.0 0.0 1.0
            0.5 0.0 1.0
            1.0 0.0 1.0
            0.0 0.3333333333333333 1.0
            0.5 0.3333333333333333 1.0
            1.0 0.3333333333333333 1.0
            0.0 0.6666666666666666 1.0
            0.5 0.6666666666666666 1.0
            1.0 0.6666666666666666 1.0
            0.0 1.0 1.0
            0.5 1.0 1.0
            1.0 1.0 1.0
        ]

        ff = FullFactorial([4, 2])
        @test doe_samples(ff, 2) == [
            0.0 0.0
            0.3333333333333333 0.0
            0.6666666666666666 0.0
            1.0 0.0
            0.0 1.0
            0.3333333333333333 1.0
            0.6666666666666666 1.0
            1.0 1.0
        ]
    end

    @testset "FractionalFactorial" begin
        frac = FractionalFactorial(["a", "b", "c", "ac"])
        @test frac.columns == ["a", "b", "c", "ac"]
        @test doe_samples(frac) == [
            0.0 0.0 0.0 1.0
            1.0 0.0 0.0 0.0
            0.0 1.0 0.0 1.0
            1.0 1.0 0.0 0.0
            0.0 0.0 1.0 0.0
            1.0 0.0 1.0 1.0
            0.0 1.0 1.0 0.0
            1.0 1.0 1.0 1.0
        ]

        frac = FractionalFactorial(["a", "b", "ba"])
        @test doe_samples(frac) == [
            0.0 0.0 1.0
            1.0 0.0 0.0
            0.0 1.0 0.0
            1.0 1.0 1.0
        ]
    end

    @testset "BoxBehnken" begin
        bb = BoxBehnken()

        @test doe_samples(bb, 2) == [0.0 0.0; 1.0 0.0; 0.0 1.0; 1.0 1.0; 0.5 0.5]
        @test doe_samples(bb, 5) == [
            0.0 0.0 0.5 0.5 0.5
            1.0 0.0 0.5 0.5 0.5
            0.0 1.0 0.5 0.5 0.5
            1.0 1.0 0.5 0.5 0.5
            0.0 0.5 0.0 0.5 0.5
            1.0 0.5 0.0 0.5 0.5
            0.0 0.5 1.0 0.5 0.5
            1.0 0.5 1.0 0.5 0.5
            0.0 0.5 0.5 0.0 0.5
            1.0 0.5 0.5 0.0 0.5
            0.0 0.5 0.5 1.0 0.5
            1.0 0.5 0.5 1.0 0.5
            0.0 0.5 0.5 0.5 0.0
            1.0 0.5 0.5 0.5 0.0
            0.0 0.5 0.5 0.5 1.0
            1.0 0.5 0.5 0.5 1.0
            0.5 0.0 0.0 0.5 0.5
            0.5 1.0 0.0 0.5 0.5
            0.5 0.0 1.0 0.5 0.5
            0.5 1.0 1.0 0.5 0.5
            0.5 0.0 0.5 0.0 0.5
            0.5 1.0 0.5 0.0 0.5
            0.5 0.0 0.5 1.0 0.5
            0.5 1.0 0.5 1.0 0.5
            0.5 0.0 0.5 0.5 0.0
            0.5 1.0 0.5 0.5 0.0
            0.5 0.0 0.5 0.5 1.0
            0.5 1.0 0.5 0.5 1.0
            0.5 0.5 0.0 0.0 0.5
            0.5 0.5 1.0 0.0 0.5
            0.5 0.5 0.0 1.0 0.5
            0.5 0.5 1.0 1.0 0.5
            0.5 0.5 0.0 0.5 0.0
            0.5 0.5 1.0 0.5 0.0
            0.5 0.5 0.0 0.5 1.0
            0.5 0.5 1.0 0.5 1.0
            0.5 0.5 0.5 0.0 0.0
            0.5 0.5 0.5 1.0 0.0
            0.5 0.5 0.5 0.0 1.0
            0.5 0.5 0.5 1.0 1.0
            0.5 0.5 0.5 0.5 0.5
            0.5 0.5 0.5 0.5 0.5
            0.5 0.5 0.5 0.5 0.5
            0.5 0.5 0.5 0.5 0.5
            0.5 0.5 0.5 0.5 0.5
            0.5 0.5 0.5 0.5 0.5
        ]

        bb = BoxBehnken(8)
        @test doe_samples(bb, 2) == [
            0.0 0.0
            1.0 0.0
            0.0 1.0
            1.0 1.0
            0.5 0.5
            0.5 0.5
            0.5 0.5
            0.5 0.5
            0.5 0.5
            0.5 0.5
            0.5 0.5
            0.5 0.5
        ]
    end

    @testset "CentralComposite" begin
        cc = CentralComposite(:face)
        @test doe_samples(cc, 3) == [
            0.0 0.0 0.0
            1.0 0.0 0.0
            0.0 1.0 0.0
            1.0 1.0 0.0
            0.0 0.0 1.0
            1.0 0.0 1.0
            0.0 1.0 1.0
            1.0 1.0 1.0
            1.0 0.5 0.5
            0.0 0.5 0.5
            0.5 1.0 0.5
            0.5 0.0 0.5
            0.5 0.5 1.0
            0.5 0.5 0.0
            0.5 0.5 0.5
        ]

        cc = CentralComposite(:inscribed)
        @test doe_samples(cc, 3) == [
            0.21132486540518708 0.21132486540518708 0.21132486540518708
            0.7886751345948129 0.21132486540518708 0.21132486540518708
            0.21132486540518708 0.7886751345948129 0.21132486540518708
            0.7886751345948129 0.7886751345948129 0.21132486540518708
            0.21132486540518708 0.21132486540518708 0.7886751345948129
            0.7886751345948129 0.21132486540518708 0.7886751345948129
            0.21132486540518708 0.7886751345948129 0.7886751345948129
            0.7886751345948129 0.7886751345948129 0.7886751345948129
            1.0 0.5 0.5
            0.0 0.5 0.5
            0.5 1.0 0.5
            0.5 0.0 0.5
            0.5 0.5 1.0
            0.5 0.5 0.0
            0.5 0.5 0.5
        ]
    end

    @testset "Variables and Generators" begin
        v, i, g, gi = UncertaintyQuantification.variables_and_generators(["a", "b", "ab"])
        @test v == "ab"
        @test i == [1, 2]
        @test g == ["ab"]
        @test gi == [3]

        v, i, g, gi = UncertaintyQuantification.variables_and_generators([
            "a", "b", "ab", "c", "ac"
        ])
        @test v == "abc"
        @test i == [1, 2, 4]
        @test g == ["ab", "ac"]
        @test gi == [3, 5]

        v, i, g, gi = UncertaintyQuantification.variables_and_generators([
            "c", "abc", "b", "a"
        ])
        @test v == "cba"
        @test i == [1, 3, 4]
        @test g == ["abc"]
        @test gi == [2]

        @test_throws ErrorException(
            "Each String in columns must hold at least one character"
        ) UncertaintyQuantification.variables_and_generators(["a", "", "b", "ab"])

        @test_throws ErrorException(
            "All combinations of columns must also be individual columns."
        ) UncertaintyQuantification.variables_and_generators(["a", "b", "ac", "ab"])
    end

    @testset "Sampling" begin
        x = RandomVariable(Uniform(-1, 1), :x)
        y = RandomVariable(Uniform(0, 1), :y)

        ff = FullFactorial([3, 3])

        samples = sample([x, y], ff)

        @test samples.x ≈ [-1, 0, 1, -1, 0, 1, -1, 0, 1] atol = 1e-12
        @test samples.y ≈ [0, 0, 0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0] atol = 1e-12

        x = RandomVariable(Normal(), :x)
        y = RandomVariable(Exponential(), :y)

        ff = FullFactorial([5, 5])

        samples = sample([x, y], ff)

        @test all(isfinite.(samples.x))
        @test all(isfinite.(samples.y))
    end
end
