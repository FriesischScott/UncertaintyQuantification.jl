using UncertaintyQuantification, DelimitedFiles

E = RandomVariable(Uniform(40, 60), :Enrich)
O = RandomVariable(Uniform(530, 690), :OuterWall)

sourcedir = joinpath(pwd() * "/../..", "demo/models")

sourcefile = "openmc_TBR.py"

numberformats = Dict(:E => ".8e")

workdir = joinpath(pwd(), "openmc_TBR")

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

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
