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

struct GibbsSampler <: AbstractBayesianMethod
    snum::Int64
    GibbsSampler(snum) = snum > 0 ? new(snum) : error("n must be greater than zero")
end

struct AdaptiveMetropolisHastings <: AbstractBayesianMethod # Adaptive Metropolis-Hastings
    snum::Integer # number of samples per chain
    cnum::Integer # number of chains
    x0::DataFrame # starting point
    C0::AbstractMatrix{<:Real} # initial proposal covariance

    function AdaptiveMetropolisHastings(
        snum::Integer, cnum::Integer, x0::DataFrame, C0::AbstractMatrix{<:Real}
    )
        errChck = (size(x0, 2) != size(C0, 1))
        errChck <<= 1
        errChck += (size(C0, 1) != size(C0, 2))
        if snum <= 0
            error("Number of samples <= 0!")
        end

        if errChck > 0
            msg = ""
            if (errChck & 1) == 1
                msg = string(msg, length(msg) > 0 ? " " : "", "Σ is not square!")
            end
            errChck >>= 1
            if (errChck & 1) == 1
                msg = string(
                    msg, length(msg) > 0 ? " " : "", "Dimension mismatch in μ and Σ!"
                )
            end
            error(msg)
        end

        return new(snum, cnum, x0, C0)
    end
end

function bayesianupdating(
    loglikelihood::Function, prior::Function, amh::AdaptiveMetropolisHastings
)
    tar(x) = exp.(loglikelihood(x)) .* prior(x)

    x0 = [copy(amh.x0) for _ in 1:(amh.cnum)]
    pmap(i -> AdaptiveMetropolisHastingsSingle!(tar, amh.snum, x0[i], amh.C0), 1:(amh.cnum))

    out = DataFrame()

    for i in 1:(amh.cnum)
        append!(out, x0[i])
    end

    return out
end

# GibbsSampler implementation
function bayesianupdating(loglikelihood::Function, prior::Function, gibbs::GibbsSampler)
    return error("Not implemented!")
end

# Adaptive Metropolis-Hastings, algorithm implemented according to Haario
# et al., "An Adaptive Metropolis Algorithm", https://doi.org/10.2307/3318737
function AdaptiveMetropolisHastingsSingle!(
    tar::Function,
    N::Integer,
    x0::DataFrame,
    C0::Matrix{<:Real},
    t0::Int64=1,
    s0::Float64=0.0001,
    burnin::Integer=0,
)
    acc = 0 # accepted Samples

    for nS in 2:(N + burnin)
        if nS < t0 || acc == 0
            C_ = C0
        else
            C_ = cov(Matrix(x0)) + s0 * C0
        end
        propD = Distributions.MvNormal(Vector(x0[nS - 1, :]), C_)
        propX = DataFrame(names(x0) .=> rand(propD))
        lastPost = tar(x0[nS - 1, :])
        propPost = tar(propX)
        #lastPost and propPost are always one-entry vectors, there should be no issue with always using the first entry?
        alpha = propPost[1] / lastPost[1]
        if alpha > 1.0 || alpha > rand(1)[1]
            append!(x0, propX)
            acc += 1
        else
            push!(x0, x0[end, :])
        end
    end
    deleteat!(x0, 1:burnin)

    return nothing
end
