# Reference: https://opensees.berkeley.edu/wiki/index.php/Simply_supported_beam_modeled_with_two_dimensional_solid_elements
using UncertaintyQuantification, DelimitedFiles, Formatting

# To run the model distributed add the desired workers and load the required packages with @everywhere
# using Distributed, Formatting
# addprocs(6; exeflags="--project")
# @everywhere using UncertaintyQuantification, DelimitedFiles

E = RandomVariable(Normal(1000, 5), :E)

# Source/Extra files are expected to be in this folder
sourcedir = joinpath(pwd(), "demo/models")

# These files will be rendere through Mustach.jl and have values injected
sourcefiles = [
    "supported-beam.tcl"
]

# These files will be copied to the working directory without injecting values
extrafiles = []

# Dictionary to map FormatSpecs (Formatting.jl) to variables
numberformats = Dict(
    :E => FormatSpec(".8e")
)

# UQ will create subfolders in here to run the solver and store the results
workdir = joinpath(pwd(), "supported-beam")

# Read output file and compute maximum (absolute) displacement
# An extractor is based the working directory for the current sample
disp = Extractor(
    base -> begin
        file = joinpath(base, "displacement.out")
        data = readdlm(file, ' ')

        return maximum(abs.(data[:, 2]))
    end,
    :disp
)

opensees = Solver(
    "/home/behrensd/src/OpenSees/bin/OpenSees",
    "", # extra args passed to the solver binary
    "supported-beam.tcl"
)

ext = ExternalModel(
    sourcedir,
    sourcefiles,
    extrafiles,
    numberformats,
    workdir,
    [disp],
    opensees
)

pf, samples = probability_of_failure(
    ext,
    df -> 0.35 .- df.disp,
    [E],
    MonteCarlo(1000)
)

println("Probability of failure: $pf")