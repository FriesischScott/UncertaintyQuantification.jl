"""
    GaussianProcess(data::DataFrame, dependendVarName::Symbol, deg::Int, dim::Int)

Creates a gaussian process prior ....

# Examples
```jldoctest
julia> data = DataFrame(x = 1:10, y = [1, 4, 10, 15, 24, 37, 50, 62, 80, 101]);

julia> rs = ResponseSurface(data, :y, 2) |> DisplayAs.withcontext(:compact => true)
ResponseSurface([0.483333, -0.238636, 1.01894], :y, [:x], 2, Monomial{Commutative{CreationOrder}, Graded{LexOrder}}[1, x₁, x₁²])
```
"""
struct GaussianProcess <: UQModel
    gp::AbstractGPs.AbstractGP
    input::Union{Vector{<:UQInput}, Vector{Symbol}}
    output::Vector{Symbol}
    inp_dim::Int
    out_dim::Int
    inp_transformer::AbstractInputTransformer # not sure if these should transform hyperparams as well
    out_transformer::AbstractOutputTransformer
end

function GaussianProcess(    
    gp::AbstractGPs.AbstractGP,
    input::Union{UQInput, Symbol},
    output::Symbol,
    inp_dim::Int,
    out_dim::Int,
    inp_transformer::AbstractInputTransformer, # not sure if these should transform hyperparams as well
    out_transformer::AbstractOutputTransformer, # leaving that for later
)
    GaussianProcess(
        gp, [input], [output], 
        inp_dim, out_dim, 
        inp_transformer, out_transformer
        )
end

function GaussianProcess(    
    gp::AbstractGPs.AbstractGP,
    input::Union{Vector{<:UQInput}, Vector{Symbol}},
    output::Symbol,
    inp_dim::Int,
    out_dim::Int,
    inp_transformer::AbstractInputTransformer, # not sure if these should transform hyperparams as well
    out_transformer::AbstractOutputTransformer, # leaving that for later
)
    GaussianProcess(
        gp, input, [output], 
        inp_dim, out_dim, 
        inp_transformer, out_transformer
        )
end

function GaussianProcess(    
    gp::AbstractGPs.AbstractGP,
    input::Union{UQInput, Symbol},
    output::Vector{Symbol},
    inp_dim::Int,
    out_dim::Int,
    inp_transformer::AbstractInputTransformer, # not sure if these should transform hyperparams as well
    out_transformer::AbstractOutputTransformer, # leaving that for later
)
    GaussianProcess(
        gp, [input], output, 
        inp_dim, out_dim, 
        inp_transformer, out_transformer
        )
end

function GaussianProcess(    
    gp::AbstractGPs.AbstractGP,
    input::Union{Vector{<:UQInput}, Vector{Symbol}},
    output::Symbol,
    inp_transformer::AbstractInputTransformer, # not sure if these should transform hyperparams as well
    out_transformer::AbstractOutputTransformer, # leaving that for later
)
    GaussianProcess(
        gp, [input], [output], 
        1, 1, 
        inp_transformer, out_transformer
        )
end

""" """
NoiseTypes = Union{
    ParameterHandling.Positive, 
    ParameterHandling.Bounded, 
    ParameterHandling.Fixed
    }

# Custom meanfunctions will break Zygote autodiff for multidimensional inputs
# Create from DataFrame
function GaussianProcess(
    data::DataFrame,
    inputs::Symbol,
    outputs::Symbol,
    build_gp::Function,
    params::NamedTuple,
    noise::NoiseTypes=positive(exp(-2.0)), # could support functions for noise as well...
    normalize_inp::Bool=false,
    normalize_out::Bool=false,
    optimizer::Union{Optim.FirstOrderOptimizer, Nothing}=nothing
)
    (inp_dim, out_dim, 
    inp_transformer, out_transformer, 
    x, y) = _handle_gp_input(
        data, inputs, outputs,
        normalize_inp, normalize_out
        )

    θ = (;
        mean_and_kernel = params,
        noise = (;noise_params = noise)
    )

    gp = build_gp_posterior(build_gp, θ, x, y, optimizer)

    return GaussianProcess(
        gp, inputs, outputs, 
        inp_dim, out_dim, 
        inp_transformer, out_transformer
        )
end

function GaussianProcess(
    data::DataFrame,
    inputs::Vector{Symbol},
    outputs::Vector{Symbol},
    build_gp::Function,
    params::NamedTuple,
    noise::NoiseTypes=positive(exp(-2.0)), # could support functions for noise as well...
    normalize_inp::Bool=false,
    normalize_out::Bool=false,
    optimizer::Union{Optim.FirstOrderOptimizer, Nothing}=nothing
)
    (inp_dim, out_dim, 
    inp_transformer, out_transformer, 
    x, y) = _handle_gp_input(
        data, inputs, outputs,
        normalize_inp, normalize_out
        )

    θ = (;
        mean_and_kernel = params,
        noise = (;noise_params = noise)
    )

    gp = build_gp_posterior(build_gp, θ, x, y, optimizer)

    return GaussianProcess(
        gp, inputs, outputs, 
        inp_dim, out_dim, 
        inp_transformer, out_transformer
        )
end

function GaussianProcess(
    data::DataFrame,
    inputs::Vector{Symbol},
    outputs::Symbol,
    build_gp::Function,
    params::NamedTuple,
    noise::NoiseTypes=positive(exp(-2.0)), # could support functions for noise as well...
    normalize_inp::Bool=false,
    normalize_out::Bool=false,
    optimizer::Union{Optim.FirstOrderOptimizer, Nothing}=nothing
)
    return GaussianProcess(
        data, inputs, [outputs],
        build_gp, params, noise, 
        normalize_inp, normalize_out, 
        optimizer
        )
end

function GaussianProcess(
    data::DataFrame,
    inputs::Symbol,
    outputs::Vector{Symbol},
    build_gp::Function,
    params::NamedTuple,
    noise::NoiseTypes=positive(exp(-2.0)), # could support functions for noise as well...
    normalize_inp::Bool=false,
    normalize_out::Bool=false,
    optimizer::Union{Optim.FirstOrderOptimizer, Nothing}=nothing
)
    return GaussianProcess(
        data, [inputs], outputs,
        build_gp, params, noise, 
        normalize_inp, normalize_out, 
        optimizer
        )
end

""" experimental design """
struct ExperimentalDesign # not sure about the name
    sim::AbstractMonteCarlo # could also allow doe
end

""" GP from uqinput and model """
# This creates a DataFrame and the calls the method above
# need to treat this differently because univariate in- and output
function GaussianProcess(
    inputs::UQInput,
    model::Vector{<:UQModel},
    outputs::Symbol,
    exp_design::ExperimentalDesign,
    build_gp::Function,
    params::NamedTuple,
    noise::NoiseTypes=positive(exp(-2.0)), # could support functions for noise as well...
    normalize_inp::Bool=false,
    normalize_out::Bool=false,
    optimizer::Union{Optim.FirstOrderOptimizer, Nothing}=nothing
)
    data = sample(inputs, exp_design.sim) # need to be able to pass experimental design
    evaluate!(model, data)

    return GaussianProcess(
        data, inputs, outputs,
        build_gp, params, noise,
        normalize_inp, normalize_out, 
        optimizer
        )
end

# need to treat this differently because univariate in- and output
function GaussianProcess(
    inputs::UQInput,
    model::UQModel,
    outputs::Symbol,
    exp_design::ExperimentalDesign,
    build_gp::Function,
    params::NamedTuple,
    noise::NoiseTypes=positive(exp(-2.0)), # could support functions for noise as well...
    normalize_inp::Bool=false,
    normalize_out::Bool=false,
    optimizer::Union{Optim.FirstOrderOptimizer, Nothing}=nothing
)
    return GaussianProcess(
        inputs, [model], outputs, exp_design,
        build_gp, params, noise,
        normalize_inp, normalize_out, 
        optimizer
        )
end

# All these cases dispatch to the same method
function GaussianProcess(
    inputs::Vector{<:UQInput},
    model::Vector{<:UQModel},
    outputs::Vector{Symbol},
    exp_design::ExperimentalDesign,
    build_gp::Function,
    params::NamedTuple,
    noise::NoiseTypes=positive(exp(-2.0)), # could support functions for noise as well...
    normalize_inp::Bool=false,
    normalize_out::Bool=false,
    optimizer::Union{Optim.FirstOrderOptimizer, Nothing}=nothing
)
    data = sample(inputs, exp_design.sim) # need to be able to pass experimental design
    evaluate!(model, data)

    return GaussianProcess(
        data, inputs, outputs,
        build_gp, params, noise,
        normalize_inp, normalize_out, 
        optimizer
        )
end

function GaussianProcess(
    inputs::Vector{<:UQInput},
    model::UQModel,
    outputs::Vector{Symbol},
    exp_design::ExperimentalDesign,
    build_gp::Function,
    params::NamedTuple,
    noise::NoiseTypes=positive(exp(-2.0)), # could support functions for noise as well...
    normalize_inp::Bool=false,
    normalize_out::Bool=false,
    optimizer::Union{Optim.FirstOrderOptimizer, Nothing}=nothing
)
    return GaussianProcess(
        inputs, [model], outputs, exp_design,
        build_gp, params, noise,
        normalize_inp, normalize_out, 
        optimizer
        )
end

function GaussianProcess(
    inputs::Vector{<:UQInput},
    model::UQModel,
    outputs::Symbol,
    exp_design::ExperimentalDesign,
    build_gp::Function,
    params::NamedTuple,
    noise::NoiseTypes=positive(exp(-2.0)), # could support functions for noise as well...
    normalize_inp::Bool=false,
    normalize_out::Bool=false,
    optimizer::Union{Optim.FirstOrderOptimizer, Nothing}=nothing
)
    return GaussianProcess(
        inputs, [model], [outputs], exp_design,
        build_gp, params, noise,
        normalize_inp, normalize_out, 
        optimizer
        )
end

function GaussianProcess(
    inputs::UQInput,
    model::Vector{<:UQModel},
    outputs::Vector{Symbol},
    exp_design::ExperimentalDesign,
    build_gp::Function,
    params::NamedTuple,
    noise::NoiseTypes=positive(exp(-2.0)), # could support functions for noise as well...
    normalize_inp::Bool=false,
    normalize_out::Bool=false,
    optimizer::Union{Optim.FirstOrderOptimizer, Nothing}=nothing
)
    return GaussianProcess(
        [inputs], model, outputs, exp_design,
        build_gp, params, noise,
        normalize_inp, normalize_out, 
        optimizer
        )
end

function GaussianProcess(
    inputs::UQInput,
    model::UQModel,
    outputs::Vector{Symbol},
    exp_design::ExperimentalDesign,
    build_gp::Function,
    params::NamedTuple,
    noise::NoiseTypes=positive(exp(-2.0)), # could support functions for noise as well...
    normalize_inp::Bool=false,
    normalize_out::Bool=false,
    optimizer::Union{Optim.FirstOrderOptimizer, Nothing}=nothing
)
    return GaussianProcess(
        [inputs], [model], outputs, exp_design,
        build_gp, params, noise,
        normalize_inp, normalize_out, 
        optimizer
        )
end

function GaussianProcess(
    inputs::Vector{<:UQInput},
    model::Vector{<:UQModel},
    outputs::Symbol,
    exp_design::ExperimentalDesign,
    build_gp::Function,
    params::NamedTuple,
    noise::NoiseTypes=positive(exp(-2.0)), # could support functions for noise as well...
    normalize_inp::Bool=false,
    normalize_out::Bool=false,
    optimizer::Union{Optim.FirstOrderOptimizer, Nothing}=nothing
)
    return GaussianProcess(
        inputs, model, [outputs], exp_design,
        build_gp, params, noise,
        normalize_inp, normalize_out, 
        optimizer
        )
end

# what should this calculate? Calculates mean for now
function evaluate!(gp::GaussianProcess, df::DataFrame) # this now gives mean and variance at input
    x = gp.inp_transformer(df)

    if gp.inp_dim != 1 || gp.out_dim != 1 # here we have to reformat the input
        x = MOInput(RowVecs(x), gp.out_dim)
    end

    y = mean(gp.gp(x))

    if gp.out_dim == 1
        y = inverse_transform(y, gp.out_transformer)
    else
        y = reshape(mean(gp.gp(x)), :, gp.out_dim)
        y = inverse_transform(y, gp.out_transformer)
    end

    insertcols!(df, (gp.output .=> eachcol(y))...)
    return nothing
end

""" maximize_logml(logml, θ, x, y, build_gp; optimizer=optimizer) """
function logml(θ, x, y, build_gp)
    gp = build_gp(ParameterHandling.value(θ.mean_and_kernel))
    f = gp(
        x, 
        only(ParameterHandling.value(θ.noise))^2 # same as in GaussianProcess...
        )
    return -logpdf(f, y)
end

function maximize_logml(logml, θ, x, y, build_gp; optimizer, maxiter=1_000)
    options = Optim.Options(; iterations=maxiter, show_trace=true)

    θ_flat, unflatten = ParameterHandling.value_flatten(θ)

    ## https://julianlsolvers.github.io/Optim.jl/stable/#user/tipsandtricks/#avoid-repeating-computations
    function fg!(F, G, θᵢ)
        if F !== nothing && G !== nothing
            val, grad = Zygote.withgradient(
                θᵢ -> logml(unflatten(θᵢ), x, y, build_gp), 
                θᵢ
                )
            G .= only(grad)
            return val
        elseif G !== nothing
            grad = Zygote.gradient(
                θᵢ -> logml(unflatten(θᵢ), x, y, build_gp), 
                θᵢ
                )
            G .= only(grad)
            return nothing
        elseif F !== nothing
            return logml(unflatten(θᵢ), x, y, build_gp)
        end
    end

    result = optimize(Optim.only_fg!(fg!), θ_flat, optimizer, options; inplace=false)

    return unflatten(result.minimizer), result
end

function build_gp_posterior(
    build_gp::Function,
    θ::NamedTuple, 
    x::AbstractArray{<:Real}, 
    y::AbstractArray{<:Real},
    optimizer::Union{Optim.FirstOrderOptimizer, Nothing}
)
    if isnothing(optimizer)
        # If no optimizer is given we just conditionalize on output
        gp = build_gp(ParameterHandling.value(θ.mean_and_kernel))
        fx = gp(x, only(ParameterHandling.value(θ.noise))^2) # this should be possible to do in a better way...
        gp = posterior(fx, y)
    else
        # Use the passed optimizer to maximize marginal log likelihood
        θ_opt, logml_ = maximize_logml(logml, θ, x, y, build_gp; optimizer=optimizer) # should I return the logml?
        gp = build_gp(ParameterHandling.value(θ_opt.mean_and_kernel))
        fx = gp(x, only(ParameterHandling.value(θ_opt.noise))^2) # this should be possible to do in a better way...
        gp = posterior(fx, y)
    end
    return gp
end
