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
        X = df_to_array(df, input_names)
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
    if !isnothing(transformer.input_transform)
        X = df_to_array(df, transformer.input_names)
        return StatsBase.transform(transformer.transform, X)
    else
        X = df_to_array(df, transformer.input_names)
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
        X = df_to_array(data, uqinput_names)
        return X
    else
        uqinput_names = names(transformer.uqinputs)
        # X is a Matrix for multiple inputs, else it is a Vector 
        X = df_to_array(df, uqinput_names)
        return X
    end
end

struct Outputtransformer <: AbstractOutputTransformer
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
        Y = df_to_array(df, output_names)
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

function (transformer::OutputTransformer)(df::DataFrame)
    if !isnothing(transformer.transform)
        Y = df_to_array(df, transformer.output_names)
        return StatsBase.transform(transformer.transform, Y)
    else
        Y = df_to_array(df, transformer.output_names)
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

function df_to_array(  # That name sucks
    df::DataFrame, 
    name::Union{Symbol, String} # do we use Strings? 
)
    return df[:, name]
end

function df_to_array(  # That name sucks
    df::DataFrame, 
    names::Union{Vector{<:Symbol}, Vector{<:String}} # do we use Strings? 
)
    # check for the case where we want a single column but the name is given in a Vector
    length(names) == 1 ? X = df_to_array(df, names[1]) : X = Matrix(df[:, names])
    return X
end