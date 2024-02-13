#===
# HPC
## Example: OpenMC TBR uncertainty

In this example, we will run OpenMC, to compute the tritium breeding ratio (TBR) uncertainty, by varying material and geometric properties. This example was taken from: https://github.com/fusion-energy/neutronics-workshop

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
