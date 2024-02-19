import UncertaintyQuantification.:run_HPC_job

sourcedir = tempdir()

# Overwrites slurm/run_HPC_job for testing
function run_HPC_job(slurm::SlurmInterface, m, path)
    dirpath = joinpath(m.workdir, path)

    p = pipeline(`bash $sourcedir/loopbash.sh 5`)
    cd(() -> run(p), dirpath)

    return nothing
end

# create source file for slurm alias
open(joinpath(sourcedir, "loopbash.sh"), "w") do input
    println(input, "#!/bin/bash")
    println(input, "")
    println(input, "total_loops=\$1")
    println(input, "")
    println(input, "for ((i=1; i<=\$total_loops; i++)); do")
    println(input, "    export SLURM_ARRAY_TASK_ID=\$i")
    println(input, "    bash slurm_array.sh")
    println(input, "done")
end

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
