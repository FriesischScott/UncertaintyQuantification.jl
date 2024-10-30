"""
Bayesian Point estimates
"""

"""
Maximum a posteriori and maximum likelihood estimates, uses optimization to find the most likely model parameters based on prior and/or likelihood.
The user can provide the method for optimization. Optim has some methods implemented, both with gradient and gradient-free.
It is possible to provide multiple initial values for optimization, MLE and MAP will return the maximum of the target function for each starting point so it is also possible to find approximations for multi-variate distributions.
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
Maximum likelihood estimate
"""
## !TODO Currently the prior is used to get information about model parameters, maybe there is a better way. In MLE the prior is not needed

struct MaximumLikelihoodBayesian <: AbstractBayesianPointEstimate

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
general function for bayesianupdating for point estimate methods. Function takes likelihood, models, estimation method and an optional prior distribution. For MLE, the prior is ignored, if no prior is given, the prior evaluation function is generated from the given prior in the struct. For both methods, the variables to be updated are inferred from the given prior in the struct, so even if the prior is not needed in MLE, it should be given as a dummy-distribution to inform the method on which variables to update. 
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

function optimize_pointestimate(optimTarget::Function, pointestimate::AbstractBayesianPointEstimate)

    method = ""

    if pointestimate.optimmethod=="LBFGS"
        method = LBFGS()
    elseif pointestimate.optimmethod=="NelderMead"
        method = NelderMead()
    else
        error("Optimization method $(pointestimate.optimmethod) is not supported in UncertaintyQuantification.jl")
    end

    if all(isinf, pointestimate.upperbounds) && all(isinf, pointestimate.lowerbounds)
        optvalues = map(x -> optimize(optimTarget, x, method), pointestimate.x0)
    else
        optvalues = map(x -> optimize(optimTarget, pointestimate.lowerbounds, pointestimate.upperbounds, x, Fminbox(method)),pointestimate.x0)
    end

end