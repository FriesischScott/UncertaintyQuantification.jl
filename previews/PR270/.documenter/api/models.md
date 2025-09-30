
# Models {#Models}

## Index {#Index}
- [`UncertaintyQuantification.Model`](#UncertaintyQuantification.Model)
- [`UncertaintyQuantification.ParallelModel`](#UncertaintyQuantification.ParallelModel)
- [`UncertaintyQuantification.UQModel`](#UncertaintyQuantification.UQModel)
- [`UncertaintyQuantification.evaluate!`](#UncertaintyQuantification.evaluate!-Tuple{ParallelModel,%20DataFrame})
- [`UncertaintyQuantification.evaluate!`](#UncertaintyQuantification.evaluate!-Tuple{Model,%20DataFrame})


## Types {#Types}
<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.UQModel' href='#UncertaintyQuantification.UQModel'><span class="jlbinding">UncertaintyQuantification.UQModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Abstract supertype for all model types


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/UncertaintyQuantification.jl#L35-L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.Model' href='#UncertaintyQuantification.Model'><span class="jlbinding">UncertaintyQuantification.Model</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Model(f::Function, name::Symbol)
```


The function `f` must accept a `DataFrame` and return the result of the model for each row in the `DataFrame` as a vector. The `name` is used to add the output to the `DataFrame`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/models/model.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.ParallelModel' href='#UncertaintyQuantification.ParallelModel'><span class="jlbinding">UncertaintyQuantification.ParallelModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ParallelModel(f::Function, name::Symbol)
```


The  `ParallelModel`  does what the `Model` does with a small difference. The function `f` is passed a `DataFrameRow` not the full `DataFrame`. If workers (through `Distributed`) are present, the rows are evaluated in parallel.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/models/model.jl#L12-L18" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Methods {#Methods}
<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.evaluate!-Tuple{Model, DataFrame}' href='#UncertaintyQuantification.evaluate!-Tuple{Model, DataFrame}'><span class="jlbinding">UncertaintyQuantification.evaluate!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
evaluate!(m::Model, df::DataFrame)
```


Calls `m.func` with `df` and adds the result to the `DataFrame` as a column `m.name`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/models/model.jl#L32-L36" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.evaluate!-Tuple{ParallelModel, DataFrame}' href='#UncertaintyQuantification.evaluate!-Tuple{ParallelModel, DataFrame}'><span class="jlbinding">UncertaintyQuantification.evaluate!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
evaluate!(m::ParallelModel, df::DataFrame)
```


Calls `m.func` for each row of `df` and adds the result to the `DataFrame` as a column `m.name`. If workers are added through `Distributed`, the rows will be evaluated in parallel.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/models/model.jl#L42-L47" target="_blank" rel="noreferrer">source</a></Badge>

</details>

