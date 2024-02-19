include("slurm_test_utils.jl")

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

    @testset "run HPC job" begin
        
        # Note, the run_HPC_job function has been overwritten in tests/slurm/slurm_test_utils.jl
        
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

        evaluate!(ext, df)
        
        @test length(readdir(readdir(ext.workdir; join=true)[1])) == 6
        @test isapprox(df.r, sqrt.(df.x .^ 2 + df.y .^ 2))        

    end
end
