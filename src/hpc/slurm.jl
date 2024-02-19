"""
	SlurmInterface(account::String, partition::String, nodes::Integer, ntasks::Integer, throttle::Integer, jobname::String, extras::Vector{String}, time::String)

When `SlurmInterface` is passed to an `ExternalModel`, model evaluations are executed using slurm job arrays.
This allows for heavier simulations or workflows to be sampled, without relying on Julia's native parallelism.
`SlurmInterface` automatically generates a slurm job array script, and Julia waits for this job to finish before extracting results.

When using `SlurmInterface`, you no longer need to load workers into Julia with `addprocs(N)`, and the requested nodes / tasks those 
required by individual model evaluations. Use `extras` to specify anything that must be preloaded for your models to be executed (for example loading modules).

The `throttle` specifies the number of simulations in the job array which are run concurrently. I.e., if you perform 
`MonteCarlo` simulation with `N=1000` samples, with `throttle=200`, it will run 1000 simulations in total, but only 200 at the same time.
Your HPC scheduler (and admin) may be unhappy if you request too many concurrent jobs.
If left empty, you scheduler's default throttle will be used.

# parameters

    account   : the account to charge the computing time to
    partition : the partition to use
    nodes     : number of nodes per job
    ntasks    : total number of cores
    throttle  : the number of jobs to be run at the same time
    jobname   : name of the job
    mempercpu : string, amount of RAM to give per cpu (default MB)
    extras    : instructions to be executed before the model is run, e.g. activating a python environment
    time      : how long each job takes, e.g. 00:10:00 for 10 minutes per sample

# Examples
```jldoctest
julia> slurm = SlurmInterface(account = "HPC_account_1", partition = "CPU_partition", nodes = 1, ntasks = 32, throttle = 200, extras = ["load python3"], time = "00:10:00")
SlurmInterface("HPC_account_1", "CPU_partition", 1, 32, 200, "UQ_array", ["load python3"], "00:10:00")


```

"""
struct SlurmInterface <: AbstractHPCScheduler
    account::String
    partition::String
    nodes::Integer
    ntasks::Integer
    throttle::Integer
    jobname::String
    mempercpu::String
    extras::Vector{String}
    time::String

    function SlurmInterface(
        account::String,
        partition::String,
        nodes::Integer,
        ntasks::Integer,
        throttle::Integer,
        jobname::String,
        mempercpu::String,
        extras::Vector{String},
        time::String,
    )
        return new(
            account, partition, nodes, ntasks, throttle, jobname, mempercpu, extras, time
        )
    end
end

function SlurmInterface(;
    account::String,
    partition::String,
    nodes::Integer,
    ntasks::Integer,
    throttle::Integer=0,
    jobname::String="UQ_array",
    mempercpu::String="",
    extras::Vector{String}=String[],
    time::String="",
)
    return SlurmInterface(
        account, partition, nodes, ntasks, throttle, jobname, mempercpu, extras, time
    )
end

function generate_slurm_job(SI, m, n, path)
    binary = m.solver.path
    source = m.solver.source
    args = m.solver.args

    extras = SI.extras

    run_command = !isempty(args) ? "$binary $args $source" : "$binary $source"
    array_command = if iszero(SI.throttle)
        "#SBATCH --array=[1-$(n)]\n"
    else
        "#SBATCH --array=[1-$(n)]%$(SI.throttle)\n"
    end

    digits = ndigits(n)

    fname = "slurm_array.sh"
    dirpath = joinpath(m.workdir, path)
    fpath = joinpath(dirpath, fname)

    open(fpath, "w") do file
        write(file, "#!/bin/bash -l\n")
        write(file, "#SBATCH -A $(SI.account)\n")
        write(file, "#SBATCH -p $(SI.partition)\n")
        write(file, "#SBATCH -J $(SI.jobname)\n")
        write(file, "#SBATCH --nodes=$(SI.nodes)\n")
        write(file, "#SBATCH --ntasks=$(SI.ntasks)\n")
        write(file, "#SBATCH --time=$(SI.time)\n")
        write(
            file, "#SBATCH --output=sample-%$(digits)a/UncertaintyQuantification-%a.out\n"
        )
        write(file, "#SBATCH --error=sample-%$(digits)a/UncertaintyQuantification-%a.err\n")
        if !isempty(SI.mempercpu)
            write(file, "#SBATCH --mem-per-cpu=$(SI.mempercpu)\n")
        end
        write(file, array_command)
        write(file, "\n\n\n")
        write(file, "#### EXTRAS ####\n")

        for extra in extras
            write(file, "$extra\n")
        end

        write(file, "\n\n\n")

        write(file, "echo ========================================================\n")
        write(file, "echo SLURM job: submitted date = `date`\n")
        write(file, "date_start=`date +%s`\n")
        write(file, "echo =========================================================\n")
        write(file, "echo Job output begins\n")
        write(file, "echo -----------------\n")
        write(file, "echo\n")
        write(file, "hostname\n")
        write(file, "echo Running with \$SLURM_NTASKS cores\n")

        write(file, "\n\n\n")

        write(file, "#### RUN COMMAND ####\n")
        write(file, "cd sample-\$(printf %0$(digits)d \$SLURM_ARRAY_TASK_ID)\n")
        write(file, "$run_command\n")

        write(file, "\n\n\n")

        write(file, "echo\n")
        write(file, "echo ---------------\n")
        write(file, "echo Job output ends\n")
        write(file, "date_end=`date +%s`\n")
        write(file, "seconds=\$((date_end-date_start))\n")
        write(file, "minutes=\$((seconds/60))\n")
        write(file, "seconds=\$((seconds-60*minutes))\n")
        write(file, "hours=\$((minutes/60))\n")
        write(file, "minutes=\$((minutes-60*hours))\n")
        write(file, "echo =========================================================\n")
        write(file, "echo SLURM job: finished date = `date`\n")
        write(
            file,
            "echo Total run time : \$hours Hours \$minutes Minutes \$seconds Seconds\n",
        )
        write(file, "echo =========================================================\n")
    end

    return nothing
end

function run_slurm_job(m, path)

    dirpath = joinpath(m.workdir, path)    

    p = pipeline(`sbatch --wait slurm_array.sh`)
    cd(() -> run(p), dirpath)

    return nothing
end