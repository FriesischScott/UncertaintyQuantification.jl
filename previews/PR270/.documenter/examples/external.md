
# External Solvers {#External-Solvers}

## OpenSees supported beam {#OpenSees-supported-beam}

In this example we will perform the reliability analysis of a supported beam using the open-source finite element software [OpenSees](https://opensees.berkeley.edu/). The example definition can be found [here](https://opensees.berkeley.edu/wiki/index.php/Simply_supported_beam_modeled_with_two_dimensional_solid_elements) We will be modeling the Young&#39;s modulus of the `ElasticIsotropic` material in the OpenSees model as a `RandomVariable`.

```julia
using UncertaintyQuantification
using DelimitedFiles # required to extract the simulation output

E = RandomVariable(Normal(1000, 5), :E) # Young's modulus
```


Source/Extra files are expected to be in this folder

```julia
sourcedir = joinpath(pwd(), "demo/models")
```


These files will be rendere through Mustache.jl and have values injected

```julia
sourcefile = "supported-beam.tcl"
```


Dictionary to map format Strings (Format.jl) to variables

```julia
numberformats = Dict(:E => ".8e")
```


UQ will create subfolders in here to run the solver and store the results

```julia
workdir = joinpath(pwd(), "supported-beam")
```


Read output file and compute maximum (absolute) displacement. The input `base` of the function used to construct the `Extractor` is the working directory for the current sample.

```julia
disp = Extractor(base -> begin
    file = joinpath(base, "displacement.out")
    data = readdlm(file, ' ')

    return maximum(abs.(data[:, 2]))
end, :disp)
```


Define the solver

```julia
opensees = Solver(
    "OpenSees", # path to OpenSees binary, here we expect it to be available on the system PATH
    "supported-beam.tcl";
    args="", # (optional) extra arguments passed to the solver
)
```


Put everything together to construct the external model

```julia
ext = ExternalModel(
    sourcedir, sourcefile, disp, opensees; workdir=workdir, formats=numberformats
)
```


Run the probability of failure analysis with 1000 samples

```julia
pf, samples = probability_of_failure(ext, df -> 0.35 .- df.disp, E, MonteCarlo(1000))

println("Probability of failure: $pf")
```


Running large simulations with an external model in series can be numerically demanding. The next example shows how to run the same example but execute the model in parallel.

## OpenSees supported beam parallel {#OpenSees-supported-beam-parallel}

load the Distributed module

```julia
using Distributed
```


add six local workers

```julia
addprocs(6)
```


The setup of the example is enclosed in the @everywhere macro. This ensures, that the main process and all worker processes have loaded the necessary modules and received all defined variables and functions. Apart from this, the setup is identical to the serial example.

```julia
@everywhere begin
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

    opensees = Solver("OpenSees", "supported-beam.tcl"; args="")

    ext = ExternalModel(
        sourcedir, sourcefile, disp, opensees; workdir=workdir, formats=numberformats
    )
end
```


The analysis is then executed on the main Julia process. When executing the model, the calls to the external solver will automatically be distributed to all available workers.

```julia
pf, samples = probability_of_failure(ext, df -> 0.35 .- df.disp, E, MonteCarlo(1000))

println("Probability of failure: $pf")
```

> 
> **Note**
> 
> Local parallelisation will quickly reach it&#39;s limits for large and complex models. Making use of  cluster thrugh the `SlurmInterface` is strongly recommended.
> 



---


_This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl)._
