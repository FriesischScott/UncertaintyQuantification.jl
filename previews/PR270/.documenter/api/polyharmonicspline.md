
# PolyharmonicSpline {#PolyharmonicSpline}

## Index {#Index}
- [`UncertaintyQuantification.PolyharmonicSpline`](#UncertaintyQuantification.PolyharmonicSpline)
- [`UncertaintyQuantification.evaluate!`](#UncertaintyQuantification.evaluate!-Tuple{PolyharmonicSpline,%20DataFrame})


## Type {#Type}
<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.PolyharmonicSpline' href='#UncertaintyQuantification.PolyharmonicSpline'><span class="jlbinding">UncertaintyQuantification.PolyharmonicSpline</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
PolyharmonicSpline(data::DataFrame, k::Int64, output::Symbol)
```


Creates a polyharmonic spline that is trained by given data.

#Examples

```julia
julia> data = DataFrame(x = 1:10, y = [1, -5, -10, -12, -8, -1, 5, 12, 23, 50]);

julia> PolyharmonicSpline(data, 2, :y) |> DisplayAs.withcontext(:compact => true)
PolyharmonicSpline([1.14733, -0.449609, 0.0140379, -1.02859, -0.219204, 0.900367, 0.00895592, 1.07145, -5.33101, 3.88628], [-112.005, 6.84443], [1.0; 2.0; … ; 9.0; 10.0;;], 2, [:x], :y)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/models/polyharmonicspline.jl#L1-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Functions {#Functions}
<details class='jldocstring custom-block' open>
<summary><a id='UncertaintyQuantification.evaluate!-Tuple{PolyharmonicSpline, DataFrame}' href='#UncertaintyQuantification.evaluate!-Tuple{PolyharmonicSpline, DataFrame}'><span class="jlbinding">UncertaintyQuantification.evaluate!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
evaluate!(ps::PolyharmonicSpline, df::DataFrame)
```


Evaluate given data using a previously contructed PolyharmonicSpline metamodel.

#Examples

```julia
julia> data = DataFrame(x = 1:10, y = [1, -5, -10, -12, -8, -1, 5, 12, 23, 50]);

julia> ps = PolyharmonicSpline(data, 2, :y);

julia> df = DataFrame( x = [2.5, 7.5, 12, 30]);

julia> evaluate!(ps, df);

julia> df.y |> DisplayAs.withcontext(:compact => true)
4-element Vector{Float64}:
  -7.75427
   8.29083
  84.4685
 260.437
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/FriesischScott/UncertaintyQuantification.jl/blob/f8bd7a9094e49042d8e9d2360393334fb1712413/src/models/polyharmonicspline.jl#L70-L92" target="_blank" rel="noreferrer">source</a></Badge>

</details>

