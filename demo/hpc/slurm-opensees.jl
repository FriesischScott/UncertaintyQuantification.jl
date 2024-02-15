using UncertaintyQuantification, DelimitedFiles


E = RandomVariable(Normal(1000, 5), :E)

sourcedir = joinpath(pwd(), "demo/models")

sourcefile = "supported-beam.tcl"

numberformats = Dict(:E => ".8e")

workdir = joinpath(pwd(), "supported-beam")

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

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
