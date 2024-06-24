struct SingleComponentMetropolisHastings <: AbstractBayesianMethod
    proposal::UnivariateDistribution
    x0::NamedTuple
    n::Int
    burnin::Int
end

function bayesianupdating(
    prior::Function,
    likelihood::Function,
    models::Vector{<:UQModel},
    mh::SingleComponentMetropolisHastings,
)
    number_of_dimensions = length(mh.x0)

    samples = DataFrame(collect(mh.x0)', collect(keys(mh.x0)))

    if !isempty(models)
        evaluate!(models, samples)
    end

    posterior = df -> likelihood(df) .* prior(df)

    rejection = 0.0

    for i in 1:(mh.n + mh.burnin - 1)
        current = DataFrame(samples[i, :])
        x = DataFrame(samples[i, :])

        for d in 1:number_of_dimensions
            x[1, d] += rand(mh.proposal)

            if !isempty(models)
                evaluate!(models, x)
            end

            α = min(1, posterior(x)[1] / posterior(current)[1])

            if α >= rand()
                current[1, d] = x[1, d]
            else
                x[1, d] = current[1, d]
                rejection += 1
            end
        end

        push!(samples, x[1, :])
    end

    rejection /= ((mh.n + mh.burnin) * number_of_dimensions)

    # discard burnin samples during return
    return samples[(mh.burnin + 1):end, :], rejection
end

function bayesianupdating(
    prior::Function,
    likelihood::Function,
    models::UQModel,
    mh::SingleComponentMetropolisHastings,
)
    return bayesianupdating(prior, likelihood, wrap(models), mh)
end

function bayesianupdating(
    prior::Function, likelihood::Function, mh::SingleComponentMetropolisHastings
)
    return bayesianupdating(prior, likelihood, UQModel[], mh)
end

struct TransitionalMarkovChainMonteCarlo <: AbstractBayesianMethod # Transitional Markov Chain Monte Carlo
    prior::Vector{RandomVariable}
    n::Int
    burnin::Int
    β::Real
    islog::Bool

    function TransitionalMarkovChainMonteCarlo(
        prior::Vector{RandomVariable}, n::Int, burnin::Int, β::Real=0.2, islog::Bool=true
    )
        if n <= 0
            error("Number of samples `n` must be positive")
        end

        return new(prior, n, burnin, β, islog)
    end
end

# TMCMC implementation
function bayesianupdating(
    prior::Function,
    likelihood::Function,
    models::Vector{<:UQModel},
    tmcmc::TransitionalMarkovChainMonteCarlo,
)
    covariance_method = LinearShrinkage(DiagonalUnitVariance(), :lw)

    rv_names = names(tmcmc.prior)

    j = 0 # iteration
    βⱼ = 0.0 # tempering

    θⱼ = sample(tmcmc.prior, tmcmc.n) # prior samples

    if !isempty(models)
        evaluate!(models, θⱼ)
    end

    S = 0.0

    while βⱼ < 1
        j += 1

        likelihood_j = tmcmc.islog ? likelihood(θⱼ) : log.(likelihood(θⱼ))

        adjust = maximum(likelihood_j)

        βⱼ⁺, wⱼ = _beta_and_weights(βⱼ, likelihood_j .- adjust)

        @debug "βⱼ" βⱼ⁺

        S += (log(mean(wⱼ)) + (βⱼ⁺ - βⱼ) * adjust)

        weights = FrequencyWeights(wⱼ ./ sum(wⱼ))

        idx = StatsBase.sample(
            collect(1:(tmcmc.n)), FrequencyWeights(wⱼ ./ sum(wⱼ)), tmcmc.n; replace=true
        )

        θⱼ⁺ = θⱼ[idx, :]

        Σⱼ = tmcmc.β^2 * cov(covariance_method, Matrix(θⱼ[:, rv_names]), weights)

        # Run inner MH algorithm

        chain = Vector{DataFrame}(undef, tmcmc.burnin + 2)

        chain[1] = copy(θⱼ⁺)

        target = if tmcmc.islog
            df -> likelihood(df) .* βⱼ⁺ .+ prior(df)
        else
            df -> log.(likelihood(df)) .* βⱼ⁺ .+ log.(prior(df))
        end

        for i in 2:(tmcmc.burnin + 2)
            next = copy(chain[i - 1])

            for (j, x) in enumerate(eachrow(next[:, rv_names]))
                next[j, rv_names] = rand(MvNormal(collect(x), Σⱼ))
            end

            # safeguard for Inf in the prior
            # !TODO: Find a cleaner way to do this
            idx_inf = findall(isinf, prior(next))

            while !isempty(idx_inf)
                for (j, x) in zip(idx_inf, (eachrow(chain[i - 1][idx_inf, rv_names])))
                    next[j, rv_names] = rand(MvNormal(collect(x), Σⱼ))
                end

                idx_inf = findall(isinf, prior(next))
            end

            if !isempty(models)
                evaluate!(models, next)
            end

            α = min.(0, target(next) .- target(chain[i - 1]))

            accept = α .>= log.(rand(length(α)))

            reject = .!accept

            next[reject, :] .= chain[i - 1][reject, :]

            chain[i] = next
        end

        θⱼ⁺ = chain[end]

        βⱼ = βⱼ⁺
        θⱼ = θⱼ⁺
    end
    return θⱼ, S
end

function bayesianupdating(
    likelihood::Function,
    models::Vector{<:UQModel},
    tmcmc::TransitionalMarkovChainMonteCarlo,
)
    prior = if tmcmc.islog
        df -> vec(
            sum(hcat(map(rv -> logpdf.(rv.dist, df[:, rv.name]), tmcmc.prior)...); dims=2),
        )
    else
        df -> vec(prod(hcat(map(rv -> pdf.(rv.dist, df[:, rv.name]), tmcmc.prior)...); dims=2))
    end

    return bayesianupdating(prior, likelihood, models, tmcmc)
end

function bayesianupdating(likelihood::Function, tmcmc::TransitionalMarkovChainMonteCarlo)
    return bayesianupdating(likelihood, UQModel[], tmcmc)
end

"""
    _beta_and_weights(β, likelihood)

Compute the next value for `β` and the nominal weights `w` using bisection.
"""
function _beta_and_weights(β::Real, L::AbstractVector{<:Real})
    low = β
    high = 2

    local x, w # Declare variables so they are visible outside the loop

    while (high - low) / middle(low, high) > 1e-6 && high > eps()
        x = middle(low, high)
        w = exp.((x - β) .* L)

        if std(w) / mean(w) > 1
            high = x
        else
            low = x
        end
    end

    if x > 1
        x = 1
        w = exp.((x - β) .* L)
    end

    return x, w
end
