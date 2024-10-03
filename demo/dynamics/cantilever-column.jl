# Reference: https://opensees.berkeley.edu/wiki/index.php?title=Time_History_Analysis_of_a_2D_Elastic_Cantilever_Column
#using UncertaintyQuantification, DelimitedFiles

# To run the model distributed add the desired workers and load the required packages with @everywhere
#using Distributed
#addprocs(6; exeflags="--project")
#@everywhere begin

using UncertaintyQuantification, DelimitedFiles

Δt = Parameter(0.02, :dt)
t = collect(0:Δt.value:10)

ω = collect(0:0.6:150)
cp = CloughPenzien(ω, 0.1, 0.8π, 0.6, 8π, 0.6)
# Ground motion
gm = SpectralRepresentation(cp, t, :gm)
gm_model = StochasticProcessModel(gm)

# Source/Extra files are expected to be in this folder
sourcedir = joinpath(pwd(), "demo/models/opensees-cantilever-column")

# These files will be rendere through Mustach.jl and have values injected
sourcefile = ["cantilever-column.tcl", "ground-motion.dat"]

# Dictionary to map format Strings (Formatting.jl) to variables
numberformats = Dict(:dt => ".8e", :gm => ".8e")

# UQ will create subfolders in here to run the solver and store the results
workdir = joinpath(pwd(), "workdir-cantilever-column")

# Read output file and compute maximum (absolute) displacement
# An extractor is based the working directory for the current sample
disp = Extractor(base -> begin
    file = joinpath(base, "displacement.out")
    data = readdlm(file, ' ')

    return maximum(abs.(data[:, 2]))
end, :disp)

opensees = Solver(
    "OpenSees", # path to OpenSees binary
    "cantilever-column.tcl";
    args="", # (optional) extra arguments passed to the solver
)

ext = ExternalModel(
    sourcedir, sourcefile, disp, opensees; workdir=workdir, formats=numberformats
)

models = [gm_model, ext]

#end #End of @everywhere

@time pf, mc_std, samples = probability_of_failure(models, df -> 250 .- df.disp, [Δt, gm], MonteCarlo(1000))

println("Probability of failure: $pf")