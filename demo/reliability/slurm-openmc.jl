# Reference: https://opensees.berkeley.edu/wiki/index.php/Simply_supported_beam_modeled_with_two_dimensional_solid_elements
using UncertaintyQuantification, DelimitedFiles

# To run the model distributed add the desired workers and load the required packages with @everywhere
# using Distributed, Formatting
# addprocs(6; exeflags="--project")
# @everywhere using UncertaintyQuantification, DelimitedFiles

E = RandomVariable(Normal(40, 60), :Enrich)
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
    name="UQ_slurm",
    account="UKAEA-AP001-CPU",
    partition="cclake",
    nodes=1,
    ntasks=1,
    batchsize=200,
    time="00:01:00",
    extras=["module load openmc", "source ~/.virtualenvs/openmc/bin/activate"],
)

ext = ExternalModel(
    sourcedir,
    sourcefile,
    TBR,
    openmc;
    workdir=workdir,
    formats=numberformats,
    slurm=slurm,
)

pf, samples = probability_of_failure(ext, df -> 1.0 .- df.TBR, [E, O], MonteCarlo(1000))

println("Probability of failure: $pf")
