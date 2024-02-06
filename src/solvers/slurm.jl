## Probably better to make a 'SlurmSolver' version of 'solver'

struct SlurmInterface

    name::String        # Does nothing
    account::String
    partition::String
    nodes::Integer
    ntasks::Integer
    batchsize::Integer
    extras::Vector{String}
    time::String

    function SlurmInterface(;name::String, 
        account::String,
        partition::String,
        nodes::Integer,
        ntasks::Integer,
        batchsize::Integer,
        extras::Vector{String},
        time::String)

        return new(name, account, partition, nodes, ntasks, batchsize, extras, time)
    end

end

function run_slurm_array(SI, m, n, path)
    #slurm_array.sh

    binary = m.solver.path
    source = m.solver.source
    args = m.solver.args

    extras = SI.extras

    run_command = !isempty(args) ? "$binary $args $source" : "$binary $source"

    digits = ndigits(n)

    fname = "slurm_array.sh"
    dirpath = joinpath(m.workdir, path)
    fpath = joinpath(dirpath, fname)
    
    open(fpath, "w") do file
        write(file, "#!/bin/bash -l\n")
        write(file, "#SBATCH -A $(SI.account)\n")
        write(file, "#SBATCH -p $(SI.partition)\n")
        write(file, "#SBATCH -J UQ_array\n")
        write(file, "#SBATCH --nodes=$(SI.nodes)\n")
        write(file, "#SBATCH --ntasks=$(SI.ntasks)\n")
        write(file, "#SBATCH --time=$(time)\n")
        write(file, "#SBATCH --output=UncertaintyQuantification_%A-%a.out\n")
        write(file, "#SBATCH --error=UncertaintyQuantification_%A-%a.err\n")
        write(file, "#SBATCH --array=[1-$(n)]%$(SI.batchsize)\n")
        write(file,"\n\n\n")
        write(file,"#### EXTRAS ####\n")

        for extra in extras
            write(file, "$extra\n")
        end

        write(file,"\n\n\n")
        write(file,"#### RUN COMMAND ####\n")
        write(file, "cd 'sample-%0$(digits)d' \$SLURM_ARRAY_TASK_ID\n")       ## Slurm array task not working
        write(file, "$run_command\n")
    end

    p = pipeline(
        Cmd(["sbatch --wait slurm_array.sh]"])
    )
    cd(() -> run(p), path)
    # cd(() -> run(p))

    return nothing
end