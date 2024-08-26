import UncertaintyQuantification.:run_slurm_array

sourcedir = tempdir()

# Overwrites slurm/run_slurm_array for testing
function run_slurm_array(
    slurm::SlurmInterface, m::ExternalModel, path::String, batch::Integer=0
)
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
