# Getting Started

Here we introduce the basic building blocks of *UncertaintyQuantification*. This includes the inputs such as `Parameter` or `RandomVariable` which will feed into any `Model` for a variety of different analyses. We will also present more advanced concepts including how to model dependencies between the inputs through copulas.

## Inputs

### Parameters

A `Parameter` is defined as a constant scalar value. In addition to value the constructor also requires a `Symbol` by which it can later be identified in the `Model`. A `Symbol` is a Julia object which is often used as a name or label. `Symbol`s are defined using the `:` prefix. Parameters represent constant deterministic values. As an example we define a `Parameter` representing the gravity of Earth.

```@example rv
using UncertaintyQuantification # hide
g = Parameter(9.81, :g)
```

Parameters are very handy when constants show up in the `Model` in multiple spaces. Instead of updating every instance in the `Model`, we can conveniently update the value by changing a single line.

### Random Variables

A `RandomVariable` is essentially a wrapper around any `UnivariateDistribution` defined in the *Distributions.jl* package [distributionsjl2021](@cite). Similarly to the `Parameter`, the second argument to the constructor is a `Symbol` acting as a unique identifier. For example, a standard gaussian random variable is defined by passing `Normal()` and `:x` as arguments.

```@example rv
x = RandomVariable(Normal(), :x)
```

A list of all possible distributions can be generated by executing `subtypes(UnivariateDistribution)` in the Julia REPL (read-eval-print loop). Note that, *Distributions* is re-exported from *UncertaintyQuantification* and no separate `using` statement is necessary. In addition, the most important methods of the `UnivariateDistribution` including `pdf`, `cdf`, and `quantile`, are also defined for the `RandomVariable`.

Random samples can be drawn from a `RandomVariable` by calling the `sample` method passing the random variable and the desired number of samples.

```@example rv
samples = sample(x, 100) # sample(x, MonteCarlo(100))
return nothing # hide
```

The `sample` method returns a `DataFrame` with the samples in a single column. When sampling from a `Vector` of random variables these invidivual columns are automatically merged into one unified `DataFrame`. By default, this will use stardard Monte Carlo simulation to obtain the samples. Alternatively, any of the quasi-Monte Carlo methods can be used instead.

```@example rv
samples = sample(x, SobolSampling(100))
samples = sample(x, LatinHypercubeSampling(100))
return nothing # hide
```

Many of the advanced simulations, e.g. line sampling or subset simulation require mappings to (and from) the standard normal space, and these are exposed through the `to_standard_normal_space!` and `to_physical_space!` methods respectively. These operate on a `DataFrame` and as such can be applied to samples directly. The transformation is done in-place, i.e. no new `DataFrame` is returned. As such, in the following example, the samples end up exactly as they were in the beginning.

```@example rv
to_standard_normal_space!(x, samples)
to_physical_space!(x, samples)
```

## Dependencies

*UncertaintyQuantification* supports modelling of dependencies through copulas. By using copulas, the modelling of the dependence structure is separated from the modelling of the univariate marginal distributions. The basis for copulas is given by Sklar's theorem [sklarFonctionsRepartitionDimensions1959](@cite). It states that any multivariate distribution $H$ in dimensions $d \geq 2$ can be separated into its marginal distributions $F_i$ and a copula function $C $.

$H(x_1,\ldots,x_2) = C(F_1(x_1),\ldots,F_d(x_d))$

For a thorough discussion of copulas, see [joeDependenceModelingCopulas2015](@cite).

In line with Sklar's theorem we build the joint distribution of two dependent random variables by separately defining the marginal distributions.

```@example copula
using UncertaintyQuantification # hide
x = RandomVariable(Normal(), :x)
y = RandomVariable(Uniform(), :y)
marginals = [x, y]
return nothing # hide
```

Next, we define the copula to model the dependence. *UncertaintyQuantification* supports Gaussian copulas for multivariate $d \geq 2$ dependence. Here, we define a Gaussian copula by passing the correlation matrix and then build the `JointDistribution` from the copula and the marginals.

```@example copula
cop = GaussianCopula([1 0.8; 0.8 1])
joint = JointDistribution(marginals, cop)
return nothing # hide
```

## Models

In this section we present the models included in *UncertaintyQuantification*. A model, in its most basic form, is a relationship between a set of input variables $x \in \mathbb{R}^{n_x}$ and an output $y \in \mathbb{R}$. Currently, most models are assumed to return single-valued outputs. However, as seen later, the `ExternalModel` is capable of extracting an arbitrary number of outputs from a single run of an external solver.

### Model

A `Model` is essentially a native Julia function operating on the previously defined inputs. Building a `Model` requires two things: a `Function`, which is internally passed a `DataFrame` containing the samples and must return a `Vector` containing the model response for each sample, and a `Symbol` which is the identifier used to add the model output into the `DataFrame`.

Suppose we wanted to define a `Model` which computes the distance from the origin of two variables $x$ and $y$ as $z$. We first define the function and then pass it to the `Model`.

```@example model
using UncertaintyQuantification, DataFrames # hide
function z(df::DataFrame)
  return @. sqrt(df.x^2 + df.y^2)
end
m = Model(z, :z)
```

An alternative for a simple model such as this, is to directly pass an anonymous function to the `Model`.

```@example model
m = Model(df -> sqrt.(df.x.^2 .+ df.y.^2), :z)
```

After defining it, a `Model` can be evaluated on a set of samples by calling the `evaluate!` method. This will add the model outcome to the `DataFrame`. Alternatively, the reponse can be obtained as a vector, by calling the `Model` as a function.

```julia
samples = sample([x, y], MonteCarlo(1000))
evaluate!(m, samples) # add to the DataFrame
output = m(samples) # return a Vector
```

However, most of the time manual evaluation of the `Model` will not be necessary as it is done internally by whichever analysis is performed.

### ParallelModel

With the basic `Model` it is up to the user to implement an efficient function which returns the model responses for all samples simultaneously. Commonly, this will involve vectorized operations as presented in the example. For more complex or longer running models, *UncertaintyQuantification* provides a simple `ParallelModel`. This model relies on the capabilites of the `Distributed` module, which is part of the standard library shipped with Julia. Without any present *workers*, the `ParallelModel` will evaluate its function in a loop for each sample. If one or more workers are present, it will automatically distribute the model evaluations. For this to work, *UncertaintyQuantification* must be loaded with the `@everywhere` macro in order to be loaded on all workers. In the following example, we first load *Distributed* and add four local workers. A simple model is then evaluated in parallel. Finally, the workers are removed.

```julia
using Distributed
addprocs(4) # add 4 local workers

@everywhere using UncertaintyQuantification

x = RandomVariable(Normal(), :x)
y = RandomVariable(Normal(), :y)

m = ParallelModel(df -> sqrt(df.x^2 .+ df.y^2), :z)

samples = sample([x, y], 1000)
evaluate!(m, samples)

rmprocs(workers()) # release the local workers
```

It is important to note, that the `ParallelModel` requires some overhead to distribute the function calls to the workers. Therefore it performs significantly slower than the standard `Model` with vectorized operations for a simple function as in this example.

By using *ClusterManagers.jl* to add the workers, the `ParallelModel` can easily be run on an existing compute cluster such as *Slurm*.

### ExternalModel

The `ExternalModel` provides interaction with almost any third-party solver. The only requirement is, that the solver uses text-based input and output files in which the values sampled from the random variables can be injected for each individual run. The output quantities are then extracated from the files generated by the solver using one (or more) `Extractor`(s). This way, the simulation techniques included in this module, can be applied to advanced models in finite element software such as *OpenSees* or *Abaqus*.

The first step in building the `ExternalModel` is to define the folder where the source files can be found as well as the working directory. Here, we assume that the source file for a simple supported beam model is located in a subdirectory of our current working directory. Similarly, the working directory for the solver is defined. In addition, we define the exact files where values need to be injected, and any extra files required. No values will be injected into the files specified as extra. In this example, no extra files are needed, so the variable is defined as an empty `String` vector.

```julia
sourcedir = joinpath(pwd(), "demo/models")
sourcefiles = ["supported-beam.tcl"]
extrafiles = String[]
workdir = joinpath(pwd(), "supported-beam")
```

Next, we must define where to inject values from the random variables and parameters into the input files. For this, we make use of the *Mustache.jl* and *Formatting.jl* modules. The values in the source file must be replaced by triple curly bracket expressions of the form `{{{ :x }}}`,  where `:x` is the identifier of the `RandomVariable` or `Parameter` to be injected. For example, to inject the Young's modulus and density of an elastic isotropic material in *OpenSees*, one could write the following.

```tcl
nDMaterial ElasticIsotropic 1 {{{ :E }}} 0.25 {{{ :rho }}}
```

This identifies where to inject the values, but not in which format. For this reason, we define a `Dict{Symbol, String}` which maps the identifiers of the inputs to a Python-style format string. In order to inject our values in scientific notation with eight digits, we use the format string `".8e"`. For any not explicitly defined `Symbol` we can include `:*` as a fallback.

```julia
formats = Dict(:E => ".8e",:rho => ".8e", :* => ".12e")
```

After formatting and injecting the values into the source file, it would look similar to this.

```tcl
nDMaterial ElasticIsotropic 1 9.99813819e+02 0.25 3.03176259e+00
```

Now that the values are injected into the source files, the next step is to extract the desired output quantities. This is done using an `Extractor`. The `Extractor` is designed similarly to the `Model` in that it takes a `Function` and a `Symbol` as its parameters. However, where a `DataFrame` is passed to the `Model`, the working directoy for the currently evaluated sample is passed to the function of the `Extractor`. The user defined function must then extract the required values from the file and return them. Here, we make use of the *DelimitedFiles* module to extract the maximum absolute displacement from the output file that *OpenSees* generated.

```julia
disp = Extractor(base -> begin
  file = joinpath(base, "displacement.out")
  data = readdlm(file, ' ')

  return maximum(abs.(data[:, 2]))
end, :disp)
```

An arbitrary number of `Extractor` functions can be defined in order to extract multiple output values from the solver.

The final step before building the model is to define the solver. The solver requires the path to the binary, and the input file. Optional command line arguments can be passed to the `Solver` through the `args` keyword. If the solver binary is not on the system path, the full path to the executable must be defined. Finally, the `ExternalModel` is assembled.

```julia
opensees = Solver(
  "OpenSees",
  "supported-beam.tcl";
  args = ""
)

ext = ExternalModel(
  sourcedir, sourcefiles, disp, opensees; formats=numberformats, workdir=workdir, extras=extrafiles
)
```

A full example of how to run a reliability analysis of a model defined in *OpenSees* can be found in the demo files.