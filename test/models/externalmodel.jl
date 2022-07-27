@testset "ExternalModel" begin
    sourcedir = tempdir()
    sourcefiles = ["in.txt"]
    extrafiles = String[]

    numberformats = Dict(:x => FormatSpec(".8e"), :* => FormatSpec(".8e"))

    workdir = map(i->(joinpath(tempdir(), "external-model-test-$(i)"), i), [true, false])

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

    ext = map(ii->ExternalModel(
        sourcedir, sourcefiles, extrafiles, numberformats, ii[1], [r], opensees, ii[2]
    ), [workdir[1], workdir[2]])

    open(joinpath(sourcedir, "in.txt"), "w") do input
        println(input, "{{{ :x }}}")
        println(input, "{{{ :y }}}")
    end

    x = RandomVariable(Uniform(0, 1), :x)
    y = RandomVariable(Uniform(0, 1), :y)

    df = sample([x, y], 1)

    for model in ext
        evaluate!(model, df)
        if model.cleanup
            @test length(readdir(readdir(model.workdir, join=true)[1])) == 0
        else
            @test length(readdir(readdir(model.workdir, join=true)[1])) != 0
        end
    end
 
    @test isapprox(df.r, sqrt.(df.x .^ 2 + df.y .^ 2))
end
