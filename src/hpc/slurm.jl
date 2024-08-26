"""
	SlurmInterface(account::String, partition::String, nodes::Integer, ntasks::Integer, throttle::Integer, jobname::String, extras::Vector{String}, time::String)

When `SlurmInterface` is passed to an `ExternalModel`, model evaluations are executed using slurm job arrays. This allows for heavier simulations or workflows to be sampled, without relying on Julia's native parallelism. `SlurmInterface` automatically generates a slurm job array script, and Julia waits for this job to finish before extracting results.

When using `SlurmInterface`, you no longer need to load workers into Julia with `addprocs(N)`, and the requested nodes / tasks those required by individual model evaluations. Use `extras` to specify anything that must be preloaded for your models to be executed (for example loading modules).

The `throttle` specifies the number of simulations in the job array which are run concurrently. I.e., if you perform
`MonteCarlo` simulation with `N=1000` samples, with `throttle=200`, it will run 1000 simulations in total, but only 200 at the same time. Your HPC scheduler (and admin) may be unhappy if you request too many concurrent jobs. If left empty, you scheduler's default throttle will be used. If your scheduler restricts the maximum size of a slurm job array, you can use `batchsize` to split the work into smaller batches.

# parameters

    account   : the account to charge the computing time to
    partition : the partition to use
    nodes     : number of nodes per job
    ntasks    : total number of cores
    throttle  : the number of jobs to be run at the same time
    batchsize : maximum size of the slurm array
    jobname   : name of the job
    mempercpu : string, amount of RAM to give per cpu (default MB)
    extras    : instructions to be executed before the model is run, e.g. activating a python environment
    time      : how long each job takes, e.g. 00:10:00 for 10 minutes per sample

# Examples
```jldoctest
julia> slurm = SlurmInterface(account = "HPC_account_1", partition = "CPU_partition", nodes = 1, ntasks = 32, throttle = 200, extras = ["load python3"], time = "00:10:00")
SlurmInterface("HPC_account_1", "CPU_partition", 1, 32, 200, "UQ_array", "", ["load python3"], "00:10:00")
```

"""
Base.@kwdef struct SlurmInterface <: AbstractHPCScheduler
    account::String
    partition::String
    nodes::Integer
    ntasks::Integer
    throttle::Integer = 0
    batchsize::Integer = 0
    jobname::String = "UQ_array"
    mempercpu::String = ""
    extras::Vector{String} = String[]
    time::String = ""
end

function setup_hpc_jobs(si::SlurmInterface, m::ExternalModel, n::Integer, datetime::String)
    if si.batchsize == 0
        setup_slurm_array(si, m::ExternalModel, n::Integer, datetime::String)
    else
        for batch in 1:ceil(Integer, n / si.batchsize)
            setup_slurm_array(si, m, n, datetime, batch)
        end
    end
end

function run_hpc_jobs(si::SlurmInterface, m::ExternalModel, n::Integer, datetime::String)
    if si.batchsize == 0
        run_slurm_array(si, m::ExternalModel, datetime::String)
    else
        for batch in 1:ceil(Integer, n / si.batchsize)
            run_slurm_array(si, m, datetime, batch)
        end
    end
end

function setup_slurm_array(
    si::SlurmInterface, m::ExternalModel, n::Integer, path::String, batch::Integer=0
)
    binary = m.solver.path
    source = m.solver.source
    args = m.solver.args

    extras = si.extras

    run_command = !isempty(args) ? "$binary $args $source" : "$binary $source"

    a, b = if batch > 0
        (batch - 1) * si.batchsize + 1, min(batch * si.batchsize, n)
    else
        1, n
    end

    array_command = if iszero(si.throttle)
        "#SBATCH --array=[$a-$b]\n"
    else
        "#SBATCH --array=[$a-$b]%$(si.throttle)\n"
    end

    digits = ndigits(n)

    fname = if iszero(batch)
        "slurm_array.sh"
    else
        "slurm_array-$batch.sh"
    end

    dirpath = joinpath(m.workdir, path)
    fpath = joinpath(dirpath, fname)

    open(fpath, "w") do file
        write(file, "#!/bin/bash -l\n")
        write(file, "#SBATCH -A $(si.account)\n")
        write(file, "#SBATCH -p $(si.partition)\n")
        write(file, "#SBATCH -J $(si.jobname)\n")
        write(file, "#SBATCH --nodes=$(si.nodes)\n")
        write(file, "#SBATCH --ntasks=$(si.ntasks)\n")
        write(
            file, "#SBATCH --output=sample-%$(digits)a/UncertaintyQuantification-%a.out\n"
        )
        write(file, "#SBATCH --error=sample-%$(digits)a/UncertaintyQuantification-%a.err\n")
        if !isempty(si.mempercpu)
            write(file, "#SBATCH --mem-per-cpu=$(si.mempercpu)\n")
        end
        if !isempty(si.time)
            write(file, "#SBATCH --time=$(si.time)\n")
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

# Slurm interface is passed to dispatch function for slurm. Perhaps there is a more elegant solution using parametric typing.
function run_slurm_array(
    si::SlurmInterface, m::ExternalModel, path::String, batch::Integer=0
)
    dirpath = joinpath(m.workdir, path)

    p = if iszero(batch)
        pipeline(`sbatch --wait slurm_array.sh`)
    else
        pipeline(`sbatch --wait slurm_array-$batch.sh`)
    end

    cd(() -> run(p), dirpath)

    return nothing
end
