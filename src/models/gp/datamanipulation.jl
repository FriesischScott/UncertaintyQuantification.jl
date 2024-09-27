abstract type AbstractInputTransformer end
abstract type AbstractOutputTransformer end

struct InputTransformer <: AbstractInputTransformer
    transform::Union{ZScoreTransform, Nothing} # there is probably a better way to do this
    input_names::Union{Symbol, Vector{<:Symbol}}
    normalize::Bool
end

# Construct from DataFrame
function InputTransformer(
    df::DataFrame, 
    input_names::Union{Symbol, Vector{<:Symbol}}, 
    normalize::Bool
)
    if normalize
        X = _dataframe_to_array(df, input_names)
        normalization = fit(ZScoreTransform, X; dims=1)
        InputTransformer(
            normalization, 
            input_names, 
            normalize
            )
    else
        InputTransformer(
            nothing, 
            input_names, 
            normalize
            )
    end
end

function (transformer::InputTransformer)(df::DataFrame)
    if !isnothing(transformer.transform)
        X = _dataframe_to_array(df, transformer.input_names)
        return StatsBase.transform(transformer.transform, X)
    else
        X = _dataframe_to_array(df, transformer.input_names)
        return X
    end
end

struct UQInputTransformer <: AbstractInputTransformer
    uqinputs::Union{UQInput, Vector{<:UQInput}}
    normalize::Bool
end

function (transformer::UQInputTransformer)(df::DataFrame)
    if transformer.normalize
        uqinput_names = names(transformer.uqinputs)
        data = df[:, uqinput_names]
        to_standard_normal_space!(transformer.uqinputs, data)
        # X is a Matrix for multiple inputs, else it is a Vector 
        X = _dataframe_to_array(data, uqinput_names)
        return X
    else
        uqinput_names = names(transformer.uqinputs)
        # X is a Matrix for multiple inputs, else it is a Vector 
        X = _dataframe_to_array(df, uqinput_names)
        return X
    end
end

struct OutputTransformer <: AbstractOutputTransformer
    transform::Union{ZScoreTransform, Nothing} # there is probably a better way to do this
    output_names::Union{Symbol, Vector{<:Symbol}}
    normalize::Bool
end

# Construct from DataFrame
function OutputTransformer(
    df::DataFrame, 
    output_names::Union{Symbol, Vector{<:Symbol}}, 
    normalize::Bool
)
    if normalize
        Y = _dataframe_to_array(df, output_names)
        normalization = fit(ZScoreTransform, Y; dims=1)
        OutputTransformer(
            normalization, 
            output_names, 
            normalize
            )
    else
        OutputTransformer(
            nothing, 
            output_names, 
            normalize
            )
    end
end

function (transformer::AbstractOutputTransformer)(df::DataFrame)
    if !isnothing(transformer.transform)
        Y = _dataframe_to_array(df, transformer.output_names)
        return StatsBase.transform(transformer.transform, Y)
    else
        Y = _dataframe_to_array(df, transformer.output_names)
        return Y
    end
end

function inverse_transform(Y::Array, transformer::AbstractOutputTransformer)
    if !isnothing(transformer.transform)
        return StatsBase.reconstruct(transformer.transform, Y)
    else
        return Y
    end
end

function _dataframe_to_array(  # That name sucks
    df::DataFrame, 
    name::Union{Symbol, String} # do we use Strings? 
)
    return df[:, name]
end

function _dataframe_to_array(  # That name sucks
    df::DataFrame, 
    names::Union{Vector{<:Symbol}, Vector{<:String}} # do we use Strings? 
)
    # check for the case where we want a single column but the name is given in a Vector
    length(names) == 1 ? X = _dataframe_to_array(df, only(names)) : X = Matrix(df[:, names])
    return X
end

function _handle_gp_input(
    data::DataFrame,
    input::Union{Symbol, UQInput},
    output::Symbol,
    normalize_inp::Bool=false,
    normalize_out::Bool=false
)
    inp_dim = 1
    out_dim = 1

    if isa(input, Symbol)
        inp_transformer = InputTransformer(data, input, normalize_inp)
    else
        inp_transformer = UQInputTransformer(input, normalize_inp)
    end
    out_transformer = OutputTransformer(data, output, normalize_out)

    # Turn DataFrame samples into X and Y arrays for GP
    x = inp_transformer(data)
    y = out_transformer(data)

    return (
        inp_dim, out_dim, 
        inp_transformer, out_transformer, 
        x, y)
end

function _handle_gp_input(
    data::DataFrame,
    inputs::Union{Vector{Symbol}, Vector{UQInput}},
    outputs::Vector{Symbol},
    normalize_inp::Bool=false,
    normalize_out::Bool=false
)
    inp_dim = length(inputs)
    out_dim = length(outputs)

    if isa(input, Symbol)
        inp_transformer = InputTransformer(data, input, normalize_inp)
    else
        inp_transformer = UQInputTransformer(input, normalize_inp)
    end
    out_transformer = OutputTransformer(data, output, normalize_out)

    # Turn DataFrame samples into X and Y arrays for GP
    X = inp_transformer(data)
    Y = out_transformer(data)
    x, y = prepare_isotopic_multi_output_data(RowVecs(X), RowVecs(Y))

    return (
        inp_dim, out_dim, 
        inp_transformer, out_transformer, 
        x, y)
end

function _handle_gp_input(
    data::DataFrame,
    inputs::Union{Vector{Symbol}, Vector{UQInput}},
    output::Symbol,
    normalize_inp::Bool=false,
    normalize_out::Bool=false
)
    return _handle_gp_input(
        data, inputs, [output],
        normalize_inp, normalize_out
        )
end

function _handle_gp_input(
    data::DataFrame,
    input::Union{Symbol, UQInput},
    outputs::Vector{Symbol},
    normalize_inp::Bool=false,
    normalize_out::Bool=false
)
    return _handle_gp_input(
        data, [input], outputs,
        normalize_inp, normalize_out
        )
end