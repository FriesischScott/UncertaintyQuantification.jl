# Reference: https://opensees.berkeley.edu/wiki/index.php/Simply_supported_beam_modeled_with_two_dimensional_solid_elements
using UncertaintyQuantification, DelimitedFiles


E = RandomVariable(Normal(1000, 5), :E)

# Source/Extra files are expected to be in this folder
sourcedir = joinpath(pathof(UncertaintyQuantification), "demo/models")

# These files will be rendere through Mustach.jl and have values injected
sourcefile = "supported-beam.tcl"

# Dictionary to map format Strings (Formatting.jl) to variables
numberformats = Dict(:E => ".8e")

# UQ will create subfolders in here to run the solver and store the results
workdir = joinpath(pwd(), "supported-beam")

# Read output file and compute maximum (absolute) displacement
# An extractor is based the working directory for the current sample
disp = Extractor(base -> begin
    file = joinpath(base, "displacement.out")
    data = readdlm(file, ' ')

    return maximum(abs.(data[:, 2]))
end, :disp)

opensees = Solver(
    "OpenSees", # path to OpenSees binary
    "supported-beam.tcl";
    args="", # (optional) extra arguments passed to the solver
)

slurm = SlurmInterface(
    jobname="UQ_slurm",
    account="EXAMPLE-0001-CPU",
    partition="cpu_parition",
    nodes=1,
    ntasks=1,
    throttle=50,
    time="00:10:00",
    extras=["module load opensees"],
)

ext = ExternalModel(
    sourcedir,
    sourcefile,
    disp,
    opensees;
    workdir=workdir,
    formats=numberformats,
    slurm=slurm,
)

pf, samples = probability_of_failure(ext, df -> 0.35 .- df.disp, E, MonteCarlo(1000))

println("Probability of failure: $pf")
