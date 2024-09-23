To run tests, from top directory of UQ.jl either:

using the package manager:

```julia
julia --project

using Pkg
Pkg.test()
```

or using the package REPL

```julia
julia --project

] test
```

or (one liner)

```
julia --project -e 'using Pkg; Pkg.test()'
```