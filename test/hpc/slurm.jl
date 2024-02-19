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

    # Function to check if (exact) line exits in file
    function isline(file, string_check)

        for (i, line) in enumerate(eachline(file))
            if (line == string_check)
                return true
            end
        end

        return false
    end

    # Checks the pattern doesn't exist anywhere
    function isnotline(file, string_check)

        for (i, line) in enumerate(eachline(file))
            if (m = match(Regex(string_check), line); m !== nothing)
                return false
            end
        end

        return true
    end
    
    @testset "Generate job" begin

        workdir = tempname()
        mkdir(workdir)

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
            scheduler = slurm
        )
        
        # run(`alias sbatch="bash $sourcedir/loopbash.sh 5"`)
        UncertaintyQuantification.generate_HPC_job(slurm, ext, size(df, 1), workdir)

        generated_file = joinpath(workdir, "slurm_array.sh")

        @test isfile(joinpath(workdir, "slurm_array.sh"))
        @test isline(generated_file, "#SBATCH -A HPC_account_1")
        @test isline(generated_file, "#SBATCH -p CPU_partition")
        @test isline(generated_file, "#SBATCH --nodes=1")
        @test isline(generated_file, "#SBATCH --ntasks=32")
        @test isline(generated_file, "#SBATCH --array=[1-5]")
        @test isnotline(generated_file, "#SBATCH --time=")
        @test isnotline(generated_file, "#SBATCH --mem-per-cpu=")

        @test isline(generated_file, "$(solver.path) $(solver.source)")

        slurm = SlurmInterface(
            jobname = "my_test_job",
            account = "HPC_account_1", 
            partition = "CPU_partition",
            time = "10:00:00",
            nodes = 2, 
            ntasks = 50,
            throttle=10,
            mempercpu="100",
            extras = ["load something", "load something else"]
        )

        UncertaintyQuantification.generate_HPC_job(slurm, ext, 100, workdir)

        @test isline(generated_file, "#SBATCH -J my_test_job")
        @test isline(generated_file, "#SBATCH --array=[1-100]%10")
        @test isline(generated_file, "#SBATCH --nodes=2")
        @test isline(generated_file, "#SBATCH --ntasks=50")
        @test isline(generated_file, "#SBATCH --time=10:00:00")
        @test isline(generated_file, "#SBATCH --mem-per-cpu=100")
        @test isline(generated_file, "load something")
        @test isline(generated_file, "load something else")


    end

    # @testset "Cleanup" begin
    #     ext = ExternalModel(
    #         sourcedir,
    #         sourcefile,
    #         radius,
    #         solver;
    #         formats=numberformats,
    #         workdir=tempname(),
    #         cleanup=true,
    #         slurm = slurm
    #     )

    #     # run(`alias sbatch="bash $sourcedir/loopbash.sh 5"`)
    #     evaluate!(ext, df)
    #     @test length(readdir(readdir(ext.workdir; join=true)[1])) == 0
    #     @test isapprox(df.r, sqrt.(df.x .^ 2 + df.y .^ 2))
    # end

    # @testset "Reuse model output" begin
    #     workdir = tempname()

    #     ext1 = ExternalModel(
    #         sourcedir, sourcefile, radius, solver; formats=numberformats, workdir=workdir, slurm = slurm
    #     )

    #     squared = Extractor(
    #         base -> begin
    #             return parse(Float64, readline(joinpath(base, "out-squared.txt")))
    #         end,
    #         :r2,
    #     )

    #     ext2 = ExternalModel(sourcedir, ["squared.jl"], squared, solver2; workdir=workdir, cleanup=true, slurm = slurm)
        
    #     # run(`alias sbatch="bash $sourcedir/loopbash.sh 5"`)
    #     evaluate!([ext1, ext2], df)
    #     @test df.r2 ≈ df.r .^ 2
    #     @test length(readdir(readdir(ext1.workdir; join=true)[1])) == 0
    #     @test length(readdir(readdir(ext2.workdir; join=true)[1])) == 0
    # end
end
