struct Parameterized{T}
    object::T
end

function (p::Parameterized)(θ)
    return apply_parameters(p.object, ParameterHandling.value(θ))
end

# """
#     parameterize(object) -> model, θ

# Turn `object` into a callable parameterized version of itself and a parameter `θ`.
# After assigning `model, θ = parameterize(object)`, calling `model(θ)` will yield the same
# `object` back. 
# """
parameterize(object) = Parameterized(object), extract_parameters(object)

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

extract_parameters(k::ConstantKernel) = ParameterHandling.positive(only(k.c))
apply_parameters(k::ConstantKernel, θ) = ConstantKernel(; c=θ)

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
# !WARNING: Incomplete
extract_parameters(t::ScaleTransform) = ParameterHandling.positive(only(t.s))
apply_parameters(::ScaleTransform, θ) = ScaleTransform(θ)

extract_parameters(t::ARDTransform) = ParameterHandling.positive(t.v)
apply_parameters(::ARDTransform, θ) = ARDTransform(θ)

# GPs
extract_parameters(f::GP) = (extract_parameters(f.mean), extract_parameters(f.kernel))
apply_parameters(f::GP, θ) = GP(
    apply_parameters(f.mean, θ[1]), apply_parameters(f.kernel, θ[2])
    )