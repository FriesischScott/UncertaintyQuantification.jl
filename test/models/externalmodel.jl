@testset "ExternalModel" begin
    sourcedir = tempdir()
    sourcefile = ["in.txt"]

    numberformats = Dict(:x => ".8e", :* => ".8e")

    radius = Extractor(
        base -> begin
            map(x -> parse(Float64, x), readlines(joinpath(base, "out.txt")))[1]
        end, :r
    )

    if Sys.iswindows()
        solver = Solver(joinpath(pwd(), "solvers/bin/radius.exe"), "", "in.txt")
        solver2 = Solver(joinpath(pwd(), "solvers/bin/squared.exe"), "", "out.txt")
    elseif Sys.isapple()
        solver = Solver(joinpath(pwd(), "solvers/bin/radius-mac"), "", "in.txt")
        solver2 = Solver(joinpath(pwd(), "solvers/bin/squared-mac"), "", "out.txt")
    else
        solver = Solver(joinpath(pwd(), "solvers/bin/radius"), "", "in.txt")
        solver2 = Solver(joinpath(pwd(), "solvers/bin/squared"), "", "out.txt")
    end

    open(joinpath(sourcedir, "in.txt"), "w") do input
        println(input, "{{{ :x }}}")
        println(input, "{{{ :y }}}")
    end

    open(joinpath(sourcedir, "extra.txt"), "w") do input
        println(input, "This is an extra file")
    end

    x = RandomVariable(Uniform(0, 1), :x)
    y = RandomVariable(Uniform(0, 1), :y)

    df = sample([x, y], 5)

    @testset "No Cleanup" begin
        ext = ExternalModel(
            sourcedir,
            sourcefile,
            radius,
            solver;
            workdir=tempname(),
            formats=numberformats,
            extras="extra.txt",
        )

        evaluate!(ext, df)
        @test length(readdir(readdir(ext.workdir; join=true)[1])) == 5
        @test "extra.txt" in
            readdir(readdir(readdir(ext.workdir; join=true)[1]; join=true)[1])
        @test isapprox(df.r, sqrt.(df.x .^ 2 .+ df.y .^ 2))
    end

    @testset "Cleanup" begin
        ext = ExternalModel(
            sourcedir,
            sourcefile,
            radius,
            solver;
            formats=numberformats,
            workdir=tempname(),
            cleanup=true,
        )
        evaluate!(ext, df)
        @test length(readdir(readdir(ext.workdir; join=true)[1])) == 0
        @test isapprox(df.r, sqrt.(df.x .^ 2 + df.y .^ 2))
    end

    @testset "Reuse model output" begin
        workdir = tempname()

        ext1 = ExternalModel(
            sourcedir, sourcefile, radius, solver; formats=numberformats, workdir=workdir
        )

        squared = Extractor(
            base -> begin
                map(x -> parse(Float64, x), readlines(joinpath(base, "out-squared.txt")))[1]
            end,
            :r2,
        )

        ext2 = ExternalModel("", "", squared, solver2; workdir=workdir, cleanup=true)

        evaluate!([ext1, ext2], df)
        @test df.r2 ≈ df.r .^ 2
        @test length(readdir(readdir(ext1.workdir; join=true)[1])) == 0
        @test length(readdir(readdir(ext2.workdir; join=true)[1])) == 0
    end
end
