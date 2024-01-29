using UncertaintyQuantification
using LinearAlgebra

function model2dof(θ::AbstractArray)

    # Initial Parameters
    # Mass
    m1 = 16531
    m2 = 16131
    M = [m1 0; 0 m2]

    # Stiffness
    kinit = 29700000
    k1 = θ[1] .* kinit
    k2 = θ[2] .* kinit
    K = [[k1 + k2 -k2]; [-k2 k2]]

    # Frequency
    D = eigvals(K, M)

    @show θ

    return sort(sqrt.(D) / 2 / pi)
end

function Likelihood(θ::AbstractArray, model::Function, Yexp::AbstractMatrix)
    sigma = 1 / 16
    Ysim = model(θ)
    J_theta = ((Ysim[1]^2 / Yexp[1]^2) - 1)^2 + ((Ysim[2]^2 / Yexp[2]^2) - 1)^2
    return exp(-J_theta / (2 * sigma^2))
end

### PRIOR
Yexp = [4.31 9.83]

function prior(x)
    return pdf(LogNormal(distribution_parameters(1.3, 1, LogNormal)...), x[1]) *
           pdf(LogNormal(distribution_parameters(0.8, 1, LogNormal)...), x[2])
end
Like(x) = Likelihood(x, model2dof, Yexp)

prop_sample() = Normal()
x0 = DataFrame(:x => 3, :y => 8)
n = 1000
burnin = 100

mh = SingleComponentMetropolisHastings(proposal, x0, n, burnin)

mh_samples = bayesianupdating(prior, Like, mh)
