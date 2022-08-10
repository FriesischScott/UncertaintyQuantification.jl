@testset "ExternalModel" begin
    sourcedir = tempdir()
    sourcefiles = ["in.txt"]
    extrafiles = String[]

    numberformats = Dict(:x => FormatSpec(".8e"), :* => FormatSpec(".8e"))

    r = Extractor(
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

    opensees = Solver(binary, "", "in.txt")

    open(joinpath(sourcedir, "in.txt"), "w") do input
        println(input, "{{{ :x }}}")
        println(input, "{{{ :y }}}")
    end

    x = RandomVariable(Uniform(0, 1), :x)
    y = RandomVariable(Uniform(0, 1), :y)

    df = sample([x, y], 1)

    @testset "No Cleanup" begin
        workdir_nocleanup = tempdir()
        ext = ExternalModel(
            sourcedir,
            sourcefiles,
            extrafiles,
            numberformats,
            workdir_nocleanup,
            [r],
            opensees,
            false,
        )
        evaluate!(ext, df)
        @test length(readdir(readdir(ext.workdir; join=true)[1])) != 0
    end

    @testset "Cleanup" begin
        workdir_cleanup = tempname()
        ext = ExternalModel(
            sourcedir,
            sourcefiles,
            extrafiles,
            numberformats,
            workdir_cleanup,
            [r],
            opensees,
            true,
        )
        evaluate!(ext, df)
        @test length(readdir(readdir(ext.workdir; join=true)[1])) == 0
    end

    @test isapprox(df.r, sqrt.(df.x .^ 2 + df.y .^ 2))
end
