@testset "ExternalModel" begin
    sourcedir = tempdir()
    sourcefiles = ["in.txt"]
    extrafiles = String[]

    numberformats = Dict(:x => ".8e", :* => ".8e")

    radius = Extractor(
        base -> begin
            map(x -> parse(Float64, x), readlines(joinpath(base, "out.txt")))[1]
        end, :r
    )

    binary = ""
    if Sys.iswindows()
        binary = joinpath(pwd(), "solvers/bin/radius.exe")
    elseif Sys.isapple()
        binary = joinpath(pwd(), "solvers/bin/radius-mac")
    else
        binary = joinpath(pwd(), "solvers/bin/radius")
    end

    solver = Solver(binary, "", "in.txt")

    open(joinpath(sourcedir, "in.txt"), "w") do input
        println(input, "{{{ :x }}}")
        println(input, "{{{ :y }}}")
    end

    x = RandomVariable(Uniform(0, 1), :x)
    y = RandomVariable(Uniform(0, 1), :y)

    df = sample([x, y], 5)

    @testset "No Cleanup" begin
        ext = ExternalModel(
            sourcedir,
            sourcefiles,
            extrafiles,
            numberformats,
            tempname(),
            [radius],
            solver,
            false,
        )
        evaluate!(ext, df)
        @test length(readdir(readdir(ext.workdir; join=true)[1])) == 5
        @test isapprox(df.r, sqrt.(df.x .^ 2 .+ df.y .^ 2))
    end

    @testset "Cleanup" begin
        ext = ExternalModel(
            sourcedir,
            sourcefiles,
            extrafiles,
            numberformats,
            tempname(),
            [radius],
            solver,
            true,
        )
        evaluate!(ext, df)
        @test length(readdir(readdir(ext.workdir; join=true)[1])) == 0
        @test isapprox(df.r, sqrt.(df.x .^ 2 + df.y .^ 2))
    end

    @testset "Reuse model output" begin
        binary2 = ""
        if Sys.iswindows()
            binary2 = joinpath(pwd(), "solvers/bin/squared.exe")
        elseif Sys.isapple()
            binary2 = joinpath(pwd(), "solvers/bin/squared-mac")
        else
            binary2 = joinpath(pwd(), "solvers/bin/squared")
        end

        solver2 = Solver(binary2, "", "out.txt")

        workdir = tempname()

        ext1 = ExternalModel(
            sourcedir,
            sourcefiles,
            extrafiles,
            numberformats,
            workdir,
            [radius],
            solver,
            false,
        )

        squared = Extractor(
            base -> begin
                map(x -> parse(Float64, x), readlines(joinpath(base, "out-squared.txt")))[1]
            end,
            :r2,
        )

        ext2 = ExternalModel(
            "",
            String[],
            extrafiles,
            Dict{Symbol,String}(),
            workdir,
            [squared],
            solver2,
            true,
        )

        evaluate!([ext1, ext2], df)
        @test df.r2 ≈ df.x .^ 2 + df.y .^ 2
    end
end
