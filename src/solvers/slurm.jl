## Probably better to make a 'SlurmSolver' version of 'solver'

struct SlurmInterface

    name::String        # Does nothing
    account::String
    nodes::Integer
    ntasks::Integer
    time::String
    partition::String
    extras::Vector{String}
    batchsize::Integer
end

function run_slurm_array(SI, m, n)
    #slurm_array.sh

    binary = ExternalModel.solver.path
    source = ExternalModel.solver.source
    args = ExternalModel.solver.args

    run_command = !isempty(args) ? "$binary $args $source" : "$binary $source"

    digits = ndigits(n)

    fname = "slurm_array.sh"
    dirpath = m.workdir
    fpath = joinpath(dirpath, fname)
    
    open(fpath, "w") do file
        write(file, "#!/bin/bash -l")
        write(file, "#SBATCH -A $(SI.account)")
        write(file, "#SBATCH -p $(SI.partition)")
        write(file, "#SBATCH -J UQ_array")
        write(file, "#SBATCH --nodes=$(SI.nodes)")
        write(file, "#SBATCH --ntasks=$(SI.ntasks)")
        write(file, "#SBATCH --time=$(time)")
        write(file, "#SBATCH --output=UncertaintyQuantification_%A-%a.out")
        write(file, "#SBATCH --error=UncertaintyQuantification_%A-%a.err")
        write(file, "#SBATCH --array=[1-$(n)]%$(SI.batch)")
        write(file,"")
        write(file,"")
        write(file,"")
        write(file,"")
        write(file,"#### EXTRAS ####")

        for extra in extras
            write(file, extra)
        end

        write(file,"#### RUN COMMAND ####")
        write(file, "cd 'sample-%0$(digits)d' $SLURM_ARRAY_TASK_ID")
        write(file, "$run_command")
    end

    p = pipeline(
        `sbatch --wait --array=[1-$(n)]%$(batch) slurm_array.sh`
    )
    cd(() -> run(p), folder)

    return nothing
end