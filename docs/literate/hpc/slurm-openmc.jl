#===
# HPC
## Slurm job arrays

When sampling large simulation models, or complicated workflows, Julia's inbuilt parallelism is insufficient. Job arrays are a useful feature of the slurm scheduler which allow you to run many similar jobs, which differ by an index (for example a sample number). This allows `UncertaintyQuantification.jl` to run heavier simulations (for example, simulations requiring multiple nodes), by offloading model sampling to a HPC machine using slurm. This way, `UncertaintyQuantification.jl` can be run on a single worker, and the HPC machine handles the rest.

For more information on job arrays, see: https://slurm.schedmd.com/job_array.html

===#

#===
## SlurmInterface

When `SlurmInterface` is passed to an `ExternalModel`, a slurm job array script is automatically generated and executed. Julia waits for this job to finish before extracting results and proceeding.

Your machine information must be passed to `SlurmInterface` for the array script to be correctly generated. For example:

===#
using UncertaintyQuantification

slurm = SlurmInterface(;
    account="HPC_account_1",
    partition="CPU_partition",
    nodes=1,
    ntasks=32,
    batchsize=50,
    extras=["load python3"],
    time="00:10:00",
)

#===

## Example: OpenMC TBR uncertainty

In this example, we will run OpenMC, to compute the tritium breeding ratio (TBR) uncertainty, by varying material and geometric properties. This example was taken from: https://github.com/fusion-energy/neutronics-workshop

===#

#md using Plots #hide
#md a = range(-5, 5; length=1000)   #hide
#md b = range(5, -5; length=1000)   #hide
#md himmelblau_f(x1, x2) = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2 #hide
#md s1 = surface(a, b, himmelblau_f; plot_title="Himmelblau's function")   #hide

#===
At first we need to create an array of random variables, that will be used when evaluating the points that our desgin produces.
It will also define the range of the function we want the design to fit.
This is also a good time to declare the function that we are working with.
===#

using UncertaintyQuantification, DelimitedFiles

E = RandomVariable(Uniform(40, 60), :Enrich)
O = RandomVariable(Uniform(530, 690), :OuterWall)

# Source/Extra files are expected to be in this folder
sourcedir = joinpath(pwd() * "/../..", "demo/models")

# These files will be rendere through Mustach.jl and have values injected
sourcefile = "openmc_TBR.py"

# Dictionary to map format Strings (Formatting.jl) to variables
numberformats = Dict(:E => ".8e")

# UQ will create subfolders in here to run the solver and store the results
workdir = joinpath(pwd(), "openmc_TBR")

# Read output file and compute maximum (absolute) displacement
# An extractor is based the working directory for the current sample
TBR = Extractor(base -> begin
    file = joinpath(base, "openmc.out")
    data = readdlm(file, ' ')

    return data[1]
end, :TBR)

openmc = Solver(
    "python3", # path to OpenSees binary
    "openmc_TBR.py";
    args="", # (optional) extra arguments passed to the solver
)

slurm = SlurmInterface(;
    jobname="UQ_slurm",
    account="EXAMPLE-0001-CPU",
    partition="cpu_parition",
    nodes=1,
    ntasks=1,
    batchsize=50,
    time="00:01:00",
    extras=["module load openmc", "source ~/.virtualenvs/openmc/bin/activate"],
)

ext = ExternalModel(
    sourcedir, sourcefile, TBR, openmc; workdir=workdir, formats=numberformats, slurm=slurm
)

pf, samples = probability_of_failure(ext, df -> 1.0 .- df.TBR, [E, O], MonteCarlo(1000))

println("Probability of failure: $pf")
