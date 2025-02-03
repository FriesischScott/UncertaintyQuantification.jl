"""
    MaximumAPosterioriBayesian(prior, optimmethod, x0; islog, lowerbounds, upperbounds)

Passed to [`bayesianupdating`](@ref) to estimate one or more maxima of the posterior distribution starting from `x0`. The optimization uses the method specified in `optimmethod`. Will calculate one estimation per point in x0. The flag `islog` specifies whether the prior and likelihood functions passed to the  [`bayesianupdating`](@ref) method are already  given as logarithms. `lowerbounds` and `upperbounds` specify optimization intervals.

Alternative constructors

```julia
    MaximumAPosterioriBayesian(prior, optimmethod, x0; islog) # `lowerbounds` = [-Inf], # `upperbounds` = [Inf]
    MaximumAPosterioriBayesian(prior, optimmethod, x0)  # `islog` = true
```
"""
struct MaximumAPosterioriBayesian <: AbstractBayesianPointEstimate

    prior::Vector{RandomVariable}
    optimmethod::String
    x0::Vector{Vector{Float64}}
    islog::Bool
    lowerbounds::Vector{Float64}
    upperbounds::Vector{Float64}

    function MaximumAPosterioriBayesian(
        prior::Vector{RandomVariable},
        optimmethod::String,
        x0::Vector{Float64};
        islog::Bool=true,
        lowerbounds::Vector{Float64} = [-Inf],
        upperbounds::Vector{Float64} = [Inf]
    )
        return MaximumAPosterioriBayesian(prior, optimmethod, [x0], islog=islog, lowerbounds=lowerbounds, upperbounds=upperbounds)
    end

    function MaximumAPosterioriBayesian(
        prior::Vector{RandomVariable},
        optimmethod::String,
        x0::Vector{Vector{Float64}};
        islog::Bool=true,
        lowerbounds::Vector{Float64} = [-Inf],
        upperbounds::Vector{Float64} = [Inf]
    )
        return new(prior, optimmethod, x0, islog, lowerbounds, upperbounds)
    end
end

"""
    MaximumLikelihoodBayesian(prior, optimmethod, x0; islog, lowerbounds, upperbounds)

Passed to [`bayesianupdating`](@ref) to estimate one or more maxima of the likelihood starting from `x0`. The optimization uses the method specified in `optimmethod`. Will calculate one estimation per point in x0. The flag `islog` specifies whether the prior and likelihood functions passed to the  [`bayesianupdating`](@ref) method are already  given as logarithms. `lowerbounds` and `upperbounds` specify optimization intervals.

Alternative constructors

```julia
    MaximumLikelihoodBayesian(prior, optimmethod, x0; islog) # `lowerbounds` = [-Inf], # `upperbounds` = [Inf]
    MaximumLikelihoodBayesian(prior, optimmethod, x0)  # `islog` = true
```
### Notes
The method uses `prior` only as information on which parameters are supposed to be optimized. The prior itself does not influence the result of the maximum likelihood estimate and can be given as a dummy distribution. For example, if two parameters `a` and `b` are supposed to be optimized, the prior could look like this

```julia
    prior = RandomVariable.(Uniform(0,1), [:a, :b])
```
"""
struct MaximumLikelihoodBayesian <: AbstractBayesianPointEstimate

    ## !TODO Currently the prior is used to get information about model parameters, maybe there is a better way. In MLE the prior is not needed

    prior::Vector{RandomVariable}
    optimmethod::String
    x0::Vector{Vector{Float64}}
    islog::Bool
    lowerbounds::Vector{Float64}
    upperbounds::Vector{Float64}

    function MaximumLikelihoodBayesian(
        prior::Vector{RandomVariable},
        optimmethod::String,
        x0::Vector{Float64};
        islog::Bool=true,
        lowerbounds::Vector{Float64} = [-Inf],
        upperbounds::Vector{Float64} = [Inf]
    )
        return MaximumLikelihoodBayesian(prior, optimmethod, [x0], islog=islog, lowerbounds=lowerbounds, upperbounds=upperbounds)
    end

    function MaximumLikelihoodBayesian(
        prior::Vector{RandomVariable},
        optimmethod::String,
        x0::Vector{Vector{Float64}};
        islog::Bool=true,
        lowerbounds::Vector{Float64} = [-Inf],
        upperbounds::Vector{Float64} = [Inf]
    )
        return new(prior, optimmethod, x0, islog, lowerbounds, upperbounds)
    end

end

"""
    bayesianupdating(likelihood, models, pointestimate; prior)

Perform bayesian updating using the given `likelihood`, `models`  and any point estimation method [`AbstractBayesianPointEstimate`](@ref).

### Notes

Method can be called with an empty Vector of models, i.e.

    bayesianupdating(likelihood, [], pointestimate)

If `prior` is not given, the method will construct a prior distribution from the prior specified in `AbstractBayesianPointEstimate.prior`.

`likelihood` is a Julia function which must be defined in terms of a `DataFrame` of samples, and must evaluate the likelihood for each row of the `DataFrame`

For example, a loglikelihood based on normal distribution using 'Data':

```julia
likelihood(df) = [sum(logpdf.(Normal.(df_i.x, 1), Data)) for df_i in eachrow(df)]
```

If a model evaluation is required to evaluate the likelihood, a vector of `UQModel`s must be passed to `bayesianupdating`. For example if the variable `x` above is the output of a numerical model.

"""
function bayesianupdating(likelihood::Function, models::Vector{<:UQModel}, pointestimate::AbstractBayesianPointEstimate; prior::Union{Function,Nothing} = nothing)

    optimTarget = setupoptimizationproblem(prior, likelihood, models, pointestimate)
    result = optimize_pointestimate(optimTarget, pointestimate)

    x = vcat(map(x -> push!(x.minimizer,-x.minimum),result))

    names = collect(p.name for p in pointestimate.prior)
    names = push!(names, :maxval)

    df = DataFrame([name => Float64[] for name in names])

    foreach(row -> push!(df,row), x)

    return df

end

# function to set up the optimization problem for MAP estimate and to generate a prior function if none is given.
function setupoptimizationproblem(prior::Union{Function, Nothing}, likelihood::Function, models::Vector{<:UQModel}, mapestimate::MaximumAPosterioriBayesian)

    # if no prior is given, generate prior function from mapestimate.prior
    if isnothing(prior)
        prior = if mapestimate.islog
            df -> vec(
                sum(hcat(map(rv -> logpdf.(rv.dist, df[:, rv.name]), mapestimate.prior)...); dims=2),
            )
        else
            df -> vec(prod(hcat(map(rv -> pdf.(rv.dist, df[:, rv.name]), mapestimate.prior)...); dims=2))
        end
    end

    target = if mapestimate.islog
        df -> prior(df) .+ likelihood(df)
    else
        df -> log.(prior(df)) .+ log.(likelihood(df))
    end

    optimTarget = x -> begin 
        names = collect(p.name for p in mapestimate.prior)
        input = DataFrame(x',names)
    
        if !isempty(models)
            evaluate!(models,input)
        end
        target(input)[1]*-1
    end

    return optimTarget

end

# function to generate the optimization problem for MLE. Note that the prior is not used, but for reasons of multiple dispatch it needs to be included
function setupoptimizationproblem(prior::Union{Function,Nothing}, likelihood::Function, models::Vector{<:UQModel}, mlestimate::MaximumLikelihoodBayesian)

    target = if mlestimate.islog
        df -> likelihood(df)
    else
        df -> log.(likelihood(df))
    end

    optimTarget = x -> begin 
        names = collect(p.name for p in mlestimate.prior)
        input = DataFrame(x',names)
    
        if !isempty(models)
            evaluate!(models,input)
        end
        target(input)[1]*-1
    end

    return optimTarget

end

# actual optimization procedure based on the point estimation method and given parameters
function optimize_pointestimate(optimTarget::Function, pointestimate::AbstractBayesianPointEstimate)

    method = LBFGS()

    if pointestimate.optimmethod=="LBFGS"
        method = LBFGS()
    elseif pointestimate.optimmethod=="NelderMead"
        method = NelderMead()
    else
        error("Optimization method $(pointestimate.optimmethod) is not supported in UncertaintyQuantification.jl. Currently supported methods are 'LBFGS' and 'NelderMead'.")
    end

    if all(isinf, pointestimate.upperbounds) && all(isinf, pointestimate.lowerbounds)
        optvalues = map(x -> optimize(optimTarget, x, method), pointestimate.x0)
    else
        optvalues = map(x -> optimize(optimTarget, pointestimate.lowerbounds, pointestimate.upperbounds, x, Fminbox(method)),pointestimate.x0)
    end

end