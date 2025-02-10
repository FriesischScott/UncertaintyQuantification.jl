abstract type AbstractInputTransform end
abstract type AbstractOutputTransform end

### No tranform, default
struct NoInputTransform <: AbstractInputTransform end

function (transform::NoInputTransform)(df::DataFrame, input_names::Union{Symbol, Vector{<:Symbol}})
    X = _dataframe_to_array(df, input_names)
    return X
end

### ZScore transform
struct InputTransform <: AbstractInputTransform
    transform::StatsBase.ZScoreTransform
    input_names::Union{Symbol, Vector{<:Symbol}}
end

# Construct from DataFrame
function InputTransform(
    df::DataFrame, 
    input_names::Union{Symbol, Vector{<:Symbol}}
)
    X = _dataframe_to_array(df, input_names)
    transform = fit(ZScoreTransform, X; dims=1)
    return InputTransform(transform, input_names)
end

function (transform::InputTransform)(df::DataFrame)
    X = _dataframe_to_array(df, transform.input_names)
    return StatsBase.transform(transform.transform, X)
end

### Standard normal transform
struct UQInputTransform <: AbstractInputTransform
    uqinputs::Union{UQInput, Vector{<:UQInput}}
end

function (transform::UQInputTransform)(df::DataFrame)
    df_copy = copy(df)
    to_standard_normal_space!(names(transform.uqinputs), df_copy)
    # X is a Matrix for multiple inputs, else it is a Vector 
    X = _dataframe_to_array(df_copy, uqinput_names)
    return X
end

### No Outputtransform, default
struct NoOutputTransform <: AbstractOutputTransform
    # we do not support multioutput yet
    output_names::Union{Symbol, Vector{<:Symbol}} # need to handle uqinput
end

function (transform::NoOutputTransform)(df::DataFrame)
    Y = _dataframe_to_array(df, transform.output_names)
    return Y
end

function inverse_transform(Y::Array, transform::NoOutputTransform)
    return Y
end

### ZScore output transform
struct OutputTransform <: AbstractOutputTransform
    transform::StatsBase.ZScoreTransform # there is probably a better way to do this
    output_names::Union{Symbol, Vector{<:Symbol}}
end

# Construct from DataFrame
function OutputTransform(
    df::DataFrame, 
    output_names::Union{Symbol, Vector{<:Symbol}}
)
    Y = _dataframe_to_array(df, output_names) # will fail if Y is not an array
    transform = fit(ZScoreTransform, Y; dims=1)
    return OutputTransform(transform, output_names)
end

function (transform::OutputTransform)(df::DataFrame)
    Y = _dataframe_to_array(df, transform.output_names)
    return StatsBase.transform(transform.transform, Y)
end

function inverse_transform(Y::Array, transform::OutputTransform)
    return StatsBase.reconstruct(transform.transform, Y)
end