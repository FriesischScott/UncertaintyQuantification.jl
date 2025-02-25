abstract type AbstractHyperparameterOptimization end

struct NoOptimization <: AbstractHyperparameterOptimization end

function optimize_hyperparameters(gp::Union{AbstractGPs.GP, NoisyGP}, x, y, opt::NoOptimization) #!TYPES
    return gp
end

struct MLE <: AbstractHyperparameterOptimization
    optimizer::Optim.FirstOrderOptimizer
    options::Optim.Options
end

MLE() = MLE(Optim.LBFGS(), Optim.Options(; iterations=1000, show_trace=true))

objective(f::Union{AbstractGPs.GP, NoisyGP}, x, y, mle::MLE) = -logpdf(f(x), y)

function optimize_hyperparameters(gp::Union{AbstractGPs.GP, NoisyGP}, x, y, mle::MLE) #!TYPES
    model, θ₀ = parameterize(gp)
    flatparams, unflatten = ParameterHandling.flatten(θ₀)

    ## https://julianlsolvers.github.io/Optim.jl/stable/#user/tipsandtricks/#avoid-repeating-computations
    function fg!(F, G, θ)
        if F !== nothing && G !== nothing
            val, grad = Zygote.withgradient(
                θ -> objective(model(unflatten(θ)), x, y, mle), 
                θ
                )
            G .= only(grad)
            return val
        elseif G !== nothing
            grad = Zygote.gradient(
                θ -> objective(model(unflatten(θ)), x, y, mle), 
                θ
                )
            G .= only(grad)
            return nothing
        elseif F !== nothing
            return objective(model(unflatten(θ)), x, y, mle)
        end
    end

    result = optimize(Optim.only_fg!(fg!), flatparams, mle.optimizer, mle.options; inplace=false)
    return model(unflatten(result.minimizer))
end