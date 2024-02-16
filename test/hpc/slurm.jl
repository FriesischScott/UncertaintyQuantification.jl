@testset "SlurmTest" begin
    sourcedir = tempdir()
    sourcefile = ["radius.jl"]

    # create source file for radius solver
    open(joinpath(sourcedir, "radius.jl"), "w") do input
        println(input, "x = {{{:x}}}")
        println(input, "y = {{{:y}}}")
        println(input, "z = sqrt(x^2+y^2)")
        println(input, "write(\"out.txt\", string(z))")
    end

    # create source file for squared solver
    open(joinpath(sourcedir, "squared.jl"), "w") do input
        println(input, "x = parse(Float64, readline(\"out.txt\"))")
        println(input, "y = x^2")
        println(input, "write(\"out-squared.txt\", string(y))")
    end

    # create source file for slurm alias

    open(joinpath(sourcedir, "loopbash.jl"), "w") do input
        println(input, "#!/bin/bash")
        println(input, "")
        println(input, "total_loops=\$1")
        println(input, "")
        println(input, "for ((i=1; i<=\$total_loops; i++)); do")
        println(input, "    export SLURM_ARRAY_TASK_ID=\$i")
        println(input, "    bash slurm_array.sh")
        println(input, "done")
    end

    numberformats = Dict(:x => ".8e", :* => ".8e")

    radius = Extractor(
        base -> begin
            return parse(Float64, readline(joinpath(base, "out.txt")))
        end, :r
    )

    binary = joinpath(Sys.BINDIR, "julia")

    solver = Solver(binary, "radius.jl")
    solver2 = Solver(binary, "squared.jl")

    open(joinpath(sourcedir, "extra.txt"), "w") do input
        println(input, "This is an extra file")
    end

    x = RandomVariable(Uniform(0, 1), :x)
    y = RandomVariable(Uniform(0, 1), :y)

    df = sample([x, y], 5)

    @testset "No Cleanup" begin

        slurm = SlurmInterface(
            account = "HPC_account_1", 
            partition = "CPU_partition", 
            nodes = 1, 
            ntasks = 32
        )

        ext = ExternalModel(
            sourcedir,
            sourcefile,
            radius,
            solver;
            workdir=tempname(),
            formats=numberformats,
            extras="extra.txt",
            slurm = slurm
        )

        run(`alias sbatch="bash $sourcedir/loopbash.sh 5"`)
        evaluate!(ext, df)

        @test length(readdir(readdir(ext.workdir; join=true)[1])) == 5
        @test "extra.txt" in
            readdir(readdir(readdir(ext.workdir; join=true)[1]; join=true)[1])
        @test df.r ≈ sqrt.(df.x .^ 2 .+ df.y .^ 2)
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

        run(`alias sbatch="bash $sourcedir/loopbash.sh 5"`)
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
                return parse(Float64, readline(joinpath(base, "out-squared.txt")))
            end,
            :r2,
        )

        ext2 = ExternalModel(sourcedir, ["squared.jl"], squared, solver2; workdir=workdir, cleanup=true)
        
        run(`alias sbatch="bash $sourcedir/loopbash.sh 5"`)
        evaluate!([ext1, ext2], df)
        @test df.r2 ≈ df.r .^ 2
        @test length(readdir(readdir(ext1.workdir; join=true)[1])) == 0
        @test length(readdir(readdir(ext2.workdir; join=true)[1])) == 0
    end
end
