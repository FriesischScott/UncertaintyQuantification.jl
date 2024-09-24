include("../test_utilities/read_write_utils.jl")

@testset "Slurm" begin
    sourcedir = tempdir()
    sourcefile = ["radius.jl"]

    # create source file for radius solver
    open(joinpath(sourcedir, "radius.jl"), "w") do input
        println(input, "x = {{{:x}}}")
        println(input, "y = {{{:y}}}")
        println(input, "z = sqrt(x^2+y^2)")
        println(input, "write(\"out.txt\", string(z))")
    end

    numberformats = Dict(:x => ".8e", :* => ".8e")

    radius = Extractor(
        base -> begin
            return parse(Float64, readline(joinpath(base, "out.txt")))
        end, :r
    )

    binary = joinpath(Sys.BINDIR, "julia")
    solver = Solver(binary, "radius.jl")

    open(joinpath(sourcedir, "extra.txt"), "w") do input
        println(input, "This is an extra file")
    end

    @testset "Setup jobs" begin
        workdir = tempname()
        mkdir(workdir)

        @testset "unbatched" begin
            slurm = SlurmInterface(;
                account=HPC_account, partition=HPC_partition, nodes=1, ntasks=32
            )

            ext = ExternalModel(
                sourcedir,
                sourcefile,
                radius,
                solver;
                workdir=tempname(),
                formats=numberformats,
                extras="extra.txt",
                scheduler=slurm,
            )

            UncertaintyQuantification.setup_hpc_jobs(slurm, ext, 100, workdir)

            generated_file = joinpath(workdir, "slurm_array.sh")

            @test isfile(generated_file)
            @test isline(generated_file, "#SBATCH -A $HPC_account")
            @test isline(generated_file, "#SBATCH -p $HPC_partition")
            @test isline(generated_file, "#SBATCH --nodes=1")
            @test isline(generated_file, "#SBATCH --ntasks=32")
            @test isline(generated_file, "#SBATCH --array=[1-100]")
            @test isnotanywhere(generated_file, "#SBATCH --time=")
            @test isnotanywhere(generated_file, "#SBATCH --mem-per-cpu=")

            @test isline(generated_file, "$(solver.path) $(solver.source)")
        end

        @testset "batched" begin
            slurm = SlurmInterface(;
                jobname="my_test_job",
                account=HPC_account,
                partition=HPC_partition,
                time="10:00:00",
                nodes=2,
                ntasks=50,
                throttle=2,
                batchsize=10,
                mempercpu="100",
                extras=["load something", "load something else"],
            )

            ext = ExternalModel(
                sourcedir,
                sourcefile,
                radius,
                solver;
                workdir=tempname(),
                formats=numberformats,
                extras="extra.txt",
                scheduler=slurm,
            )

            UncertaintyQuantification.setup_hpc_jobs(slurm, ext, 100, workdir)

            for batch in 1:10
                generated_file = joinpath(workdir, "slurm_array-$batch.sh")

                a = (batch - 1) * 10 + 1
                b = min(batch * 10, 100)

                @test isfile(generated_file)
                @test isline(generated_file, "#SBATCH -J my_test_job")
                @test isline(generated_file, "#SBATCH --array=[$a-$b]%2")
                @test isline(generated_file, "#SBATCH --nodes=2")
                @test isline(generated_file, "#SBATCH --ntasks=50")
                @test isline(generated_file, "#SBATCH --time=10:00:00")
                @test isline(generated_file, "#SBATCH --mem-per-cpu=100")
                @test isline(generated_file, "load something")
                @test isline(generated_file, "load something else")
            end
        end
    end

    @testset "Run jobs" begin

        # Note, the run_HPC_job function has been overwritten in tests/test_utilities/slurm_test_utils.jl

        @testset "unbatched" begin
            slurm = SlurmInterface(;
                account=HPC_account,
                partition=HPC_partition,
                nodes=1,
                ntasks=1,
                time="00:01:00",
            )

            ext = ExternalModel(
                sourcedir,
                sourcefile,
                radius,
                solver;
                workdir="hpc/test_dir",
                formats=numberformats,
                extras="extra.txt",
                scheduler=slurm,
            )

            x = RandomVariable(Uniform(0, 1), :x)
            y = RandomVariable(Uniform(0, 1), :y)

            df = sample([x, y], 10)

            evaluate!(ext, df)

            @test isapprox(df.r, sqrt.(df.x .^ 2 + df.y .^ 2))
        end

        @testset "batched" begin
            slurm = SlurmInterface(;
                account=HPC_account,
                partition=HPC_partition,
                nodes=1,
                ntasks=1,
                batchsize=4,
                time="00:01:00",
            )

            ext = ExternalModel(
                sourcedir,
                sourcefile,
                radius,
                solver;
                workdir="hpc/test_dir",
                formats=numberformats,
                extras="extra.txt",
                scheduler=slurm,
            )

            x = RandomVariable(Uniform(0, 1), :x)
            y = RandomVariable(Uniform(0, 1), :y)

            df = sample([x, y], 10)

            evaluate!(ext, df)

            @test isapprox(df.r, sqrt.(df.x .^ 2 + df.y .^ 2))
        end
    end
end
