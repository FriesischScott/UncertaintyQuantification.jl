#===
# High Performance Computing
## OpenMC TBR uncertainty

In this example, we will run OpenMC, to compute the tritium breeding ratio (TBR) uncertainty, by varying material and geometric properties. This example was taken from: https://github.com/fusion-energy/neutronics-workshop

At first we need to create an array of random variables, that will be used when evaluating the points that our desgin produces.
It will also define the range of the function we want the design to fit.
This is also a good time to declare the function that we are working with.
===#

#===
Here we will vary the models Li6 enrichment, and radius of the tungsten layer
===#
using UncertaintyQuantification, DelimitedFiles

E = RandomVariable(Uniform(0.3, 0.7), :E)
R1 = RandomVariable(Uniform(8, 14), :R1)

# Source/Extra files are expected to be in this folder
sourcedir = joinpath(pwd(), "demo/models")

# These files will be rendere through Mustach.jl and have values injected
sourcefile = "openmc_TBR.py"

# Dictionary to map format Strings (Formatting.jl) to variables
numberformats = Dict(:E => ".8e")

# UQ will create subfolders in here to run the solver and store the results
workdir = joinpath(pwd(), "openmc_TBR")

# Read output file and compute the Tritium breeding ratio
# An extractor is based the working directory for the current sample
TBR = Extractor(base -> begin
    file = joinpath(base, "openmc.out")
    data = readdlm(file, ' ')

    return data[1]
end, :TBR)

#===
In this example, an OpenMC model is built and run using the python API. We therefore specify the python3 command, and the python file to run
===#

openmc = Solver(
    "python3", # path to python3 binary
    "openmc_TBR.py";
    args="", # (optional) extra arguments passed to the solver
)

#===
Slurm interface, passing requiring machine information. Note in `extras`, we specify the commands we rquire to run the model (for example, loading modulefiles or data)
===#

slurm = SlurmInterface(;
    jobname="UQ_slurm",
    account="EXAMPLE-0001-CPU",
    partition="cpu_parition",
    nodes=1,
    ntasks=1,
    throttle=50,
    time="00:01:00",    # Time per simulation
    extras=["module load openmc", "source ~/.virtualenvs/openmc/bin/activate"],
)

# Define external model, passing slurm interface
ext = ExternalModel(
    sourcedir, sourcefile, TBR, openmc; workdir=workdir, formats=numberformats, slurm=slurm
)

# Specify a limitstate function, negative value consititutes failure. Here we are interested in P(TBR <= 1).
function limitstate(df)
    return  reduce(vcat, df.TBR) .- 1 
end


# Run a Monte Carlo simulation
@time pf, cov, samples = probability_of_failure(ext, limitstate, [E, R1], MonteCarlo(5000))
println("Probability of failure: $pf")


#jl using StatsBase

#jl TBR = samples.TBR #hide
#jl TBR_mean = mean(TBR) #hide
#jl TBR_std  = std(TBR) #hide
#jl lower_quantile = quantile(TBR, 0.025) #hide
#jl upper_quantile = quantile(TBR, 0.975) #hide
#jl println("TBR mean: $TBR_mean, TBR std: $TBR_std, TBR 95%: [$lower_quantile, $upper_quantile]") #hide

#===
Gives `pf = 0.0064`, `cov = 0.1762` for `5000` samples. Also `TBR mean = 1.2404114576444423`, `TBR std = 0.10000460056126671`, and `TBR 95%: [1.0379178446904211, 1.4130216792418262]`. 

We can also try with Subset simulation
===#

subset = SubSetInfinity(800, 0.1, 10, 0.5)

@time pf, cov, samples = probability_of_failure(ext, limitstate, [E, R1], subset)
println("Probability of failure: $pf")

#===

Giving `pf = 0.0053`, `cov = 0.2138` for `2400` samples.

===#