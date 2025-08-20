using ParameterHandling
using AbstractGPs


struct Parameterized{T}
    object::T
end

function (p::Parameterized)(θ)
    return apply_parameters(p.object, ParameterHandling.value(θ))
end

parameterize(object) = Parameterized(object), extract_parameters(object)

extract_parameters(f::GP) = (extract_parameters(f.mean), extract_parameters(f.kernel))
function apply_parameters(f::GP, θ)
    return GP(apply_parameters(f.mean, θ[1]), apply_parameters(f.kernel, θ[2]))
end

extract_parameters(m::ConstMean) = m.c
apply_parameters(::ConstMean, θ) = ConstMean(θ)

KernelsWithoutParameters = Union{SEKernel,Matern32Kernel,Matern52Kernel,WhiteKernel}

extract_parameters(::T) where {T<:KernelsWithoutParameters} = nothing
apply_parameters(k::T, θ) where {T<:KernelsWithoutParameters} = k

# ------------------------------------------

extract_parameters(k::TransformedKernel) = (extract_parameters(k.kernel), extract_parameters(k.transform))
apply_parameters(k::TransformedKernel, θ) = TransformedKernel(
    apply_parameters(k.kernel, θ[1]), apply_parameters(k.transform, θ[2])
    )

extract_parameters(k::PeriodicKernel) = ParameterHandling.positive(k.r)
apply_parameters(::PeriodicKernel, θ) = PeriodicKernel(; r=θ)

extract_parameters(t::ARDTransform) = ParameterHandling.positive(t.v)
apply_parameters(::ARDTransform, θ) = ARDTransform(θ)

extract_parameters(t::LinearTransform) = t.A
apply_parameters(::LinearTransform, θ) = LinearTransform(θ)

extract_parameters(k::KernelSum) = map(extract_parameters, k.kernels)
apply_parameters(k::KernelSum, θ) = KernelSum(map(apply_parameters, k.kernels, θ))

kernel = PeriodicKernel(2)
kernel = PeriodicKernel(1)

transform = LinearTransform(rand(2,2))

extract_parameters(kernel)
extract_parameters(transform)

A = rand(2,2)
kernel = SqExponentialKernel() ∘ LinearTransform(A)
kernel = (SqExponentialKernel() ∘ ARDTransform([1.0, 1.0])) ⊗ (SqExponentialKernel() ∘ LinearTransform(A))
gp = GP(0.0, kernel)

params = extract_parameters(gp)
model, θ = parameterize(gp)
θ_flat, unflatten = ParameterHandling.flatten(θ)
gp_model = model(unflatten(θ_flat))


struct FixedTransform{T} <: Transform
    component::T
end

struct FixedKernel{T} <: Kernel
    component::T
end

fixed(t::Transform) = FixedTransform(t)
fixed(k::Kernel) = FixedKernel(k)

extract_parameters(c::FixedTransform) = ParameterHandling.fixed(extract_parameters(c.component))
apply_parameters(c::FixedTransform, θ) = apply_parameters(c.component, θ)

extract_parameters(c::FixedKernel) = ParameterHandling.fixed(extract_parameters(c.component))
apply_parameters(c::FixedKernel, θ) = apply_parameters(c.component, θ)

A = rand(2,2)
kernel = fixed((SqExponentialKernel() ∘ LinearTransform(A)) + (SqExponentialKernel() ∘ ARDTransform([1.0, 1.0])))
gp = GP(0.0, kernel)

params = extract_parameters(gp)
model, θ = parameterize(gp)
θ_flat, unflatten = ParameterHandling.flatten(θ)
gp_model = model(unflatten(θ_flat))

function collect_concrete_kernels(T::Type)
    result = Set{Type}()

    function recurse(t)
        for s in subtypes(t)
            if isabstracttype(s)
                recurse(s)  # dive into abstract types
            else
                push!(result, s)  # collect concrete type
            end
        end
    end

    recurse(T)
    return collect(result)
end

# Retrieve all kernel types
all_kernels = collect_concrete_kernels(KernelFunctions.Kernel)
for k in all_kernels
    println(k)
end
println("Found $(length(all_kernels)) concrete kernel types:")

all_transforms = all_concrete_subtypes(KernelFunctions.Transform)
println("Found $(length(all_transforms)) concrete transformation types.")

# Retrieve all transformation types
transform_types = all_transform_types()
println("Found $(length(transform_types)) transformation types.")

for KT in all_kernel_types
    try
        obj = KT()                     # maybe you need default constructors
        extract_parameters(obj)
    catch e
        println("Missing extract_parameters for $KT: $e")
    end
end