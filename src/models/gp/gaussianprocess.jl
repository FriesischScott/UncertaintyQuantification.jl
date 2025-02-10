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
    gp::AbstractGPs.FiniteGP
    input::Vector{Symbol}
    output::Symbol
    input_transform::AbstractInputTransform # not sure if these should transform hyperparams as well
    output_transform::AbstractOutputTransform
end

""" ---------------------------------------------------------------------------------------------------- """
# names = propertynames(data[:, Not(output)])
## DATAFRAME INPUT OUTPUT STUFF
function GaussianProcess(
    gp::AbstractGPs.GP,
    data::DataFrame,
    output::Symbol,
    input_transform::NoInputTransform,
    output_transform::NoOutputTransform,
    optimization::NoOptimization
) # should we use keyword args?

    return
end

function GaussianProcess(
    gp::AbstractGPs.GP,
    data::DataFrame,
    output::Symbol,
) # should we use keyword args?
    return GaussianProcess(
        gp,
        data,
        output,
        NoInputTransform(),
        NoOutputTransform(),
        NoOptimization()
    )
end

function GaussianProcess(
    gp::AbstractGPs.GP,
    data::DataFrame,
    output::Symbol,
    input_transform::NoInputTransform,
    output_transform::NoOutputTransform,
    optimization::AbstractHyperparameterOptimization
) # should we use keyword args?

    return
end

function GaussianProcess(
    gp::AbstractGPs.GP,
    data::DataFrame,
    output::Symbol;
    input_transform::AbstractInputTransformer,
    output_transform::AbstractOutputTransformer,
    optimization::AbstractHyperparameterOptimization
)
    return
end

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

""" ---------------------------------------------------------------------------------------------------- """
## UQMODEL STUFF
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

    # random_inputs = filter(i -> isa(i, RandomUQInput), inputs)
    # random_input_names = names(random_inputs) 

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
        gp, inputs.name, outputs, 
        inp_dim, out_dim, 
        inp_transformer, out_transformer
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
