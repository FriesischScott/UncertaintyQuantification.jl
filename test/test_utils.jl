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


