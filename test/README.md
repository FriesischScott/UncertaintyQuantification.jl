To run tests, from top directory of UQ.jl either:

using the package manager:

```julia
shell> julia --project

julia> using Pkg
julia> Pkg.test()
```

or using the package REPL

```julia
shell> julia --project

julia>] test
```

or (one liner)

```
shell> julia --project -e 'using Pkg; Pkg.test()'
```

## Testing Slurm on HPC

We only test slurm on Linux distributions. For non-HPC systems, we use a dummy command that replaces sbatch in `test/test_utilities`. This may need to be added to PATH if slurm tests fail.

You may also test the Slurm interface on your HPC machine using actual slurm:

```julia
shell> julia --project

julia> using Pkg 
julia> Pkg.test(;test_args=["HPC", "YOUR_ACCOUNT", "YOUR_PARTITION"])
```

To test different configurations of the HPC interface, we launch 4 lightweight (> 1 minute) test job arrays using 1 task per job.