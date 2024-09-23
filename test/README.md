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