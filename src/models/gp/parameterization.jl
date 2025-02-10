struct Parameterized{T}
    object::T
end

function (p::Parameterized)(θ)
    return apply_parameters(p.object, ParameterHandling.value(θ))
end

"""
    parameterize(object) -> model, θ

Turn `object` into a callable parameterized version of itself and a parameter `θ`.
After assigning `model, θ = parameterize(object)`, calling `model(θ)` will yield the same
`object` back. 
"""
parameterize(object) = Parameterized(object), extract_parameters(object)

# Mean functions
extract_parameters(::ZeroMean) = nothing
apply_parameters(m::ZeroMean, θ) = m

extract_parameters(m::ConstMean) = m.c
apply_parameters(::ConstMean, θ) = ConstMean(θ)

# Simple kernels
KernelsWithoutParameters = Union{SEKernel,Matern32Kernel,Matern52Kernel,WhiteKernel}

extract_parameters(::T) where {T<:KernelsWithoutParameters} = nothing
apply_parameters(k::T, θ) where {T<:KernelsWithoutParameters} = k

extract_parameters(k::PeriodicKernel) = ParameterHandling.positive(only(k.r))
apply_parameters(::PeriodicKernel, θ) = PeriodicKernel(; r=[θ])

extract_parameters(k::RationalQuadraticKernel) = ParameterHandling.positive(only(k.α))
apply_parameters(k::RationalQuadraticKernel, θ) = RationalQuadraticKernel(; α=θ, metric=k.metric)

# Composite kernels
extract_parameters(k::KernelSum) = map(extract_parameters, k.kernels)
apply_parameters(k::KernelSum, θ) = KernelSum(map(apply_parameters, k.kernels, θ))

extract_parameters(k::KernelProduct) = map(extract_parameters, k.kernels)
apply_parameters(k::KernelProduct, θ) = KernelProduct(map(apply_parameters, k.kernels, θ))

extract_parameters(k::TransformedKernel) = (extract_parameters(k.kernel), extract_parameters(k.transform))
apply_parameters(k::TransformedKernel, θ) = TransformedKernel(
    apply_parameters(k.kernel, θ[1]), apply_parameters(k.transform, θ[2])
    )

extract_parameters(k::ScaledKernel) = (extract_parameters(k.kernel), ParameterHandling.positive(only(k.σ²)))
apply_parameters(k::ScaledKernel, θ) = ScaledKernel(apply_parameters(k.kernel, θ[1]), θ[2])

# Transforms
extract_parameters(t::ScaleTransform) = ParameterHandling.positive(only(t.s))
apply_parameters(::ScaleTransform, θ) = ScaleTransform(θ)

# # Likelihoods
# extract_parameters(::BernoulliLikelihood) = nothing
# apply_parameters(l::BernoulliLikelihood, θ) = l
# _isequal(l1::T, l2::T) where {T<:BernoulliLikelihood} = true

# extract_parameters(::PoissonLikelihood) = nothing
# apply_parameters(l::PoissonLikelihood, θ) = l
# _isequal(l1::T, l2::T) where {T<:PoissonLikelihood} = true

# GPs
extract_parameters(f::GP) = (extract_parameters(f.mean), extract_parameters(f.kernel))
apply_parameters(f::GP, θ) = GP(
    apply_parameters(f.mean, θ[1]), apply_parameters(f.kernel, θ[2])
    )

# extract_parameters(f::LatentGP) = (extract_parameters(f.f), extract_parameters(f.lik))
# function apply_parameters(f::LatentGP, θ)
#     return LatentGP(apply_parameters(f.f, θ[1]), apply_parameters(f.lik, θ[2]), f.Σy)
# end

# # Approximations
# const SVA = SparseVariationalApproximation

# function extract_parameters(sva::SVA, fixed_inducing_points::Bool)
#     fz_par = fixed_inducing_points ? nothing : collect(sva.fz.x)
#     q_par = extract_parameters(sva.q)
#     return (fz_par, q_par)
# end

# function apply_parameters(sva::SVA, θ)
#     fz = isnothing(θ[1]) ? sva.fz : sva.fz.f(θ[1])
#     q = apply_parameters(sva.q, θ[2])
#     return SVA(fz, q)
# end

# variational_gaussian(n::Int, T=Float64) = MvNormal(zeros(T, n), Matrix{T}(I, n, n))

# # Distributions
# extract_parameters(d::MvNormal) = (d.μ, ParameterHandling.positive_definite(d.Σ))
# apply_parameters(::MvNormal, θ) = MvNormal(θ[1], θ[2])
# _isequal(d1::MvNormal, d2::MvNormal) = isapprox(d1.μ, d1.μ) && isapprox(d1.Σ, d2.Σ)

# Custom wrappers
struct NoisyGP{T<:GP,Tn<:Real}
    gp::T
    obs_noise::Tn
end

(gp::NoisyGP)(x) = gp.gp(x, gp.obs_noise)

with_gaussian_noise(gp::GP, obs_noise::Real) = NoisyGP(gp, obs_noise)

extract_parameters(f::NoisyGP) = (
    extract_parameters(f.gp), 
    ParameterHandling.positive(f.obs_noise, exp, 1e-6)
    )
apply_parameters(f::NoisyGP, θ) = NoisyGP(apply_parameters(f.gp, θ[1]), θ[2])

# function _isequal(f1::NoisyGP, f2::NoisyGP)
#     return _isequal(f1.gp, f2.gp) && isapprox(f1.obs_noise, f2.obs_noise)
# end

# struct SVGP{T<:LatentGP,Ts<:SVA}
#     lgp::T
#     sva::Ts
#     fixed_inducing_points::Bool
# end

# SVGP(lgp, sva; fixed_inducing_points) = SVGP(lgp, sva, fixed_inducing_points)

# function extract_parameters(f::SVGP)
#     return (extract_parameters(f.lgp), extract_parameters(f.sva, f.fixed_inducing_points))
# end

# function apply_parameters(f::SVGP, θ)
#     lgp = apply_parameters(f.lgp, θ[1])
#     sva = apply_parameters(f.sva, θ[2])
#     return SVGP(lgp, SVA(lgp(sva.fz.x).fx, sva.q), f.fixed_inducing_points)
# end

# costfunction(svgp::SVGP, data) = -elbo(svgp.sva, svgp.lgp(data.x), data.y)