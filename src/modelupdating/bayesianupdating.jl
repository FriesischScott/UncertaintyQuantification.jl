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

struct TMCMC <: AbstractBayesianMethod # Transitional Markov Chain Monte Carlo
    sample_prior::Function
    n::Int
    burnin::Int
    β2::Real
    TMCMC(n) = n > 0 ? new(n) : error("n must be greater than zero")
end

# TMCMC implementation
function bayesianupdating(
    prior::Function, likelihood::Function, models::Vector{<:UQModel}, tmcmc::TMCMC
)
    αj = [0]
    j = 1
    θj = tmcmc.sample_prior(n)
    log_likelihood_j = zeros(n, 1)
    log_evidence = [0]

    θ_dim = size(θj, 2)

    while αj[j] < 1
        print("Start iteration $j")

        # Parallel Computation of Likelihood
        temp_likelihood_j = pmap(likelihood, eachrow(θj))
        log_likelihood_j = log.(temp_likelihood_j)
        log_likelihood_j = reduce(vcat, log_likelihood_j)

        # BisectionMethod
        low_αj = αj[j]
        up_αj = 2

        log_likelihood_adjust = maximum(log_likelihood_j)
        ωj_adjust = zeros(size(log_likelihood_j, 1), 1)
        while (up_αj - low_αj) / ((up_αj + low_αj) / 2) > 1e-6
            α_new = (up_αj + low_αj) / 2
            weight_test = exp(α_new - αj[j]) .* (log_likelihood_j .- log_likelihood_adjust)
            if (std(weight_test) / mean(weight_test)) > 1
                up_αj = α_new
            else
                low_αj = α_new
            end
            ωj_adjust = weight_test
        end
        push!(αj, minimum(1, α_new))

        # Compute plausibility weights
        ωj = ωj_adjust / sum(ωj_adjust)

        # Compute log evidence
        push!(
            log_evidence, log(mean(ωj_adjust)) + (αj[j + 1] - αj[j]) * log_likelihood_adjust
        )

        # MH porposal pdf: covariance matrix and weighted mean

    end
    return error("Not implemented!")
end
