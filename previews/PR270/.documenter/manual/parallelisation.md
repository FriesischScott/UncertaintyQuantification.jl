
# Parallelisation {#Parallelisation}

_UncertaintyQuantification_ provides several possibilities to parallelise the execution of a model. Which method to use strongly depends on the model and the available hardware. The simplest way to run a model in parallel is by using the [`ParallelModel`](/manual/parallelisation#ParallelModel).

## ParallelModel {#ParallelModel}

With the basic [`Model`](/manual/gettingstarted#Model) it is up to the user to implement an efficient function which returns the model responses for all samples simultaneously. Commonly, this will involve vectorized operations. For more complex or longer running models, _UncertaintyQuantification_ provides a simple [`ParallelModel`](/manual/parallelisation#ParallelModel). This model relies on the capabilites of the `Distributed` module, which is part of the standard library shipped with Julia. Without any present _workers_, the [`ParallelModel`](/manual/parallelisation#ParallelModel) will evaluate its function in a loop for each sample. If one or more workers are present, it will automatically distribute the model evaluations. For this to work, _UncertaintyQuantification_ must be loaded with the `@everywhere` macro in order to be loaded on all workers. In difference to the standard `Model` each call to the function used in the  `ParallelModel` is passed a `DataFrameRow` instead of the full `DataFrame`.

In the following example, we first load `Distributed` and add four local workers. A simple model is then evaluated in parallel. Finally, the workers are removed.

```julia
using Distributed
addprocs(4) # add 4 local workers

# setup of the model and inputs
@everywhere begin

using UncertaintyQuantification

x = RandomVariable(Normal(), :x)
y = RandomVariable(Normal(), :y)

m = ParallelModel(df -> sqrt(df.x^2 .+ df.y^2), :z)

end

samples = sample([x, y], 1000)
evaluate!(m, samples)

rmprocs(workers()) # release the local workers
```


It is important to note, that the [`ParallelModel`](/manual/parallelisation#ParallelModel) requires some overhead to distribute the function calls to the workers. Therefore it performs significantly slower than the standard [`Model`](/manual/gettingstarted#Model) with vectorized operations for a simple function as in this example. The method of executing a model in parallel presented above also applies to the `ExternalModel`.

::: tip Note

For heavier external models the use of parallel compute clusters through the [`SlurmInterface`](/api/slurm#UncertaintyQuantification.SlurmInterface) is recommended. See [High Performance Computing](/examples/hpc#High-Performance-Computing).

:::
