using DifferentiationInterface

abstract type AbstractHyperparameterOptimization end

struct NoOptimization <: AbstractHyperparameterOptimization end

struct MaximumLikelihoodEstimation <: AbstractHyperparameterOptimization
    optimizer::Optim.FirstOrderOptimizer
    options::Optim.Options
end

MaximumLikelihoodEstimation() = MaximumLikelihoodEstimation(
    Optim.LBFGS(),
    Optim.Options(; iterations=1000, show_trace=false)
)

function optimize_hyperparameters(
    gp::Union{AbstractGPs.GP, NoisyGP}, 
    x, 
    y, 
    opt::NoOptimization
)
    return gp
end

objective(
    f::Union{AbstractGPs.GP, NoisyGP}, 
    x, 
    y, 
    mle::MaximumLikelihoodEstimation
) = -logpdf(f(x), y)

function optimize_hyperparameters(
    gp::Union{AbstractGPs.GP, NoisyGP}, 
    x, 
    y, 
    mle::MaximumLikelihoodEstimation
) #!TYPES
    model, θ₀ = parameterize(gp)
    θ₀_flat, unflatten = ParameterHandling.flatten(θ₀)

    result = optimize(
        θ -> objective(model(unflatten(θ)), x, y, mle), 
        θ₀_flat, 
        mle.optimizer, mle.options; 
        autodiff=AutoZygote()
        )
    return model(unflatten(result.minimizer))
end