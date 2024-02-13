# HPC
## Slurm job arrays

When sampling large simulation models, or complicated workflows, Julia's inbuilt parallelism is sometimes insufficient. Job arrays are a useful feature of the slurm scheduler which allow you to run many similar jobs, which differ by an index (for example a sample number). This allows `UncertaintyQuantification.jl` to run heavier simulations (for example, simulations requiring multiple nodes), by offloading model sampling to an HPC machine using slurm. This way, `UncertaintyQuantification.jl` can be run on a single worker, and the HPC machine handles the rest.

For more information on job arrays, see: [job arrays](https://slurm.schedmd.com/job_array.html)


## SlurmInterface

When `SlurmInterface` is passed to an `ExternalModel`, a slurm job array script is automatically generated and executed. Julia waits for this job to finish before extracting results and proceeding.

Your machine information must be passed to `SlurmInterface` for the array script to be correctly generated. For example:

```@example hpc
using UncertaintyQuantification

slurm = SlurmInterface(;
    account="HPC_account_1",
    partition="CPU_partition",
    nodes=1,
    ntasks=32,
    batchsize=50,
    extras=["load python3"],
    time="01:00:00",
)
```

Here `account` is your account (provided by your HPC admin/PI), and `partition` specifies the queue that jobs will be submitted to (ask admin if unsure). `nodes` and `ntasks` are the number of nodes and cpus that your individual simulations requires. Depending on your HPC machine, each node has a specific number of cpus. If your application requires more cpus than are available per node, you can use multiple nodes.

`time` specifies the maximum time that each simulation will be run for, before being killed.

!!! warning "Individual model runs VS overall batch"
    `nodes`, `ntasks`, and `time` are parameters required for each _individual_ model evaluation, not the entire batch. For example, if you are running a large FEM simulation that requires 100 cpus to evaluate one sample, and your HPC machine has 50 cpus per node, you would specify `nodes = 2` and `ntasks = 100`.


!!! note "Compiling with MPI"
    If your model requires multiple `nodes`, it may be best to compile your application with MPI, if your model allows for it. Please check your application's documentation for compiling with MPI.

`extras` are commands that will be executed before you model is run, for example loading any module files or data your model requires. Multiple commands can be passed: `extras = ["load python", "python3 get_data.py"]`.

`batchsize` specifies the job array task throttle: the number of samples that will be run at any given time. For example, when running a `MonteCarlo` simulation with 2000 samples, and `batchsize = 50`, 2000 model evaluations will be run in total, but only 50 at the same time. If left empty, your scheduler's default throttle will be used.


<!-- ## Useage

!!! note
    See [eamples/HPC](../examples/hpc.md) for a detailed example

Let's make a simple external model, consisting of a python script. -->