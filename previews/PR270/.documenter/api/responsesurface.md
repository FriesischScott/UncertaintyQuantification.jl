
# ResponseSurface {#ResponseSurface}

## Index {#Index}


## Type {#Type}
<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.ResponseSurface' href='#UncertaintyQuantification.ResponseSurface'><span class="jlbinding">UncertaintyQuantification.ResponseSurface</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ResponseSurface(data::DataFrame, dependendVarName::Symbol, deg::Int, dim::Int)
```


Creates a response surface using polynomial least squares regression with given degree.

**Examples**

```julia
julia> data = DataFrame(x = 1:10, y = [1, 4, 10, 15, 24, 37, 50, 62, 80, 101]);

julia> rs = ResponseSurface(data, :y, 2)
ResponseSurface([0.48333333333332457, -0.23863636363636026, 1.0189393939393936], :y, [:x], 2, Monomials.Monomial[1, x1, x1²])
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/models/responsesurface.jl#L1-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Functions {#Functions}
<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.evaluate!-Tuple{ResponseSurface, DataFrame}' href='#UncertaintyQuantification.evaluate!-Tuple{ResponseSurface, DataFrame}'><span class="jlbinding">UncertaintyQuantification.evaluate!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
evaluate!(rs::ResponseSurface, data::DataFrame)
```


evaluating data by using a previously trained ResponseSurface.

**Examples**

```julia
julia> data = DataFrame(x = 1:10, y = [1, 4, 10, 15, 24, 37, 50, 62, 80, 101]);

julia> rs = ResponseSurface(data, :y, 2);

julia> df = DataFrame(x = [2.5, 11, 15]);

julia> evaluate!(rs, df);

julia> df.y |> DisplayAs.withcontext(:compact => true)
3-element Vector{Float64}:
   6.25511
 121.15
 226.165
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/models/responsesurface.jl#L39-L62" target="_blank" rel="noreferrer">source</a></Badge>

</details>

