
#===

## OpenSees supported beam parallel
===#

# load the Distributed module
using Distributed

# add six local workers
addprocs(6)

# The setup of the example is enclosed in the @everywhere macro. This ensures, that the main process and all worker processes have loaded the necessary modules and received all defined variables and functions.

@everywhere begin

using UncertaintyQuantification, DelimitedFiles

E = RandomVariable(Normal(1000, 5), :E) # Young's modulus

# Source/Extra files are expected to be in this folder
sourcedir = joinpath(pwd(), "demo/models")

# These files will be rendere through Mustache.jl and have values injected
sourcefile = "supported-beam.tcl"

# Dictionary to map format Strings (Format.jl) to variables
numberformats = Dict(:E => ".8e")

# UQ will create subfolders in here to run the solver and store the results
workdir = joinpath(pwd(), "supported-beam")

# Read output file and compute maximum (absolute) displacement
# An extractor is passed the working directory for the current sample as `base`
disp = Extractor(base -> begin
    file = joinpath(base, "displacement.out")
    data = readdlm(file, ' ')

    return maximum(abs.(data[:, 2]))
end, :disp)

# Define the solver
opensees = Solver(
    "OpenSees", # path to OpenSees binary, here we expect it to be available on the system PATH
    "supported-beam.tcl";
    args="", # (optional) extra arguments passed to the solver
)

# Put everything together to construct the external model
ext = ExternalModel(
    sourcedir, sourcefile, disp, opensees; workdir=workdir, formats=numberformats
)

end

# The analysis is then executed on the main Julia process. When executing the model, the calls to the external solver will automatically be distributed to all available workers.

# Run the probability of failure analysis with 1000 samples
pf, samples = probability_of_failure(ext, df -> 0.35 .- df.disp, E, MonteCarlo(1000))

println("Probability of failure: $pf")

#===
!!! note
    Local parallelisation will quickly reach it's limits for large and complex models. Making use of  cluster thrugh the [`SlurmInterface`](@ref) is strongly recommended.
===#
