struct StandardizeInput
    flag::Bool
end

StandardizeInput() = StandardizeInput(true)

struct StandardizeOutput
    flag::Bool
end

StandardizeOutput() = StandardizeOutput(true)

# DataTransform
abstract type AbstractInputTransform end
abstract type AbstractOutputTransform end

struct DataTransform
    input::AbstractInputTransform
    output::AbstractOutputTransform
end

function DataTransform(
    data::DataFrame,
    instandard::StandardizeInput,
    input::Union{Symbol, Vector{<:Symbol}},
    outstandard::StandardizeOutput,
    output::Symbol 
)
    intransform = instandard.flag ? InputTransform(data, input) : NoInputTransform(input)
    outtransform = outstandard.flag ? OutputTransform(data, output) : NoOutputTransform(output)
    return DataTransform(intransform, outtransform)
end

function DataTransform(
    data::DataFrame,
    instandard::StandardizeInput,
    input::Union{UQInput, Vector{<:UQInput}},
    outstandard::StandardizeOutput,
    output::Symbol 
)
    intransform = instandard.flag ? UQInputTransform(input) : NoInputTransform(input)
    outtransform = outstandard.flag ? OutputTransform(data, output) : NoOutputTransform(output)
    return DataTransform(intransform, outtransform)
end

# Utility functions to handle DataFrames
function _dataframe_to_array(
    df::DataFrame, 
    name::Symbol 
)
    return df[:, name]
end

function _dataframe_to_array(
    df::DataFrame, 
    names::Vector{<:Symbol}
)
    # check for the case where we want a single column but the name is given in a Vector
    length(names) == 1 ? X = _dataframe_to_array(df, only(names)) : X = RowVecs(Matrix(df[:, names]))
    return X
end

### No transform, default
struct NoInputTransform <: AbstractInputTransform 
    input::Union{Symbol, Vector{<:Symbol}}
end

NoInputTransform(input::Union{UQInput, Vector{<:UQInput}}) = NoInputTransform(names(input))

(transform::NoInputTransform)(df::DataFrame) = _dataframe_to_array(df, transform.input)


### ZScore transform
struct InputTransform <: AbstractInputTransform
    transform::StatsBase.ZScoreTransform
    input::Union{Symbol, Vector{<:Symbol}}
end

# Construct from DataFrame
function InputTransform(
    df::DataFrame, 
    input::Union{Symbol, Vector{<:Symbol}}
)
    X = _dataframe_to_array(df, input)
    transform = fit(ZScoreTransform, X; dims=1)
    return InputTransform(transform, input)
end

function (transform::InputTransform)(df::DataFrame)
    X = _dataframe_to_array(df, transform.input)
    return StatsBase.transform(transform.transform, X)
end

### Standard normal transform
struct UQInputTransform <: AbstractInputTransform
    uqinput::Union{UQInput, Vector{<:UQInput}}
end

function (transform::UQInputTransform)(df::DataFrame)
    df_copy = copy(df)
    uqinput_names = names(transform.uqinput)
    to_standard_normal_space!(transform.uqinput, df_copy)
    # X is a Matrix for multiple inputs, else it is a Vector 
    X = _dataframe_to_array(df_copy, uqinput_names)
    return X
end

### No Outputtransform, default
struct NoOutputTransform <: AbstractOutputTransform 
    output::Symbol
end

(transform::NoOutputTransform)(df::DataFrame) = _dataframe_to_array(df, transform.output)
inverse_transform(Y::Array, transform::NoOutputTransform) = Y

### ZScore output transform
struct OutputTransform <: AbstractOutputTransform
    transform::StatsBase.ZScoreTransform # there is probably a better way to do this
    output::Symbol
end

# Construct from DataFrame
function OutputTransform(
    df::DataFrame, 
    output::Symbol
)
    Y = _dataframe_to_array(df, output) # will fail if Y is not an array
    transform = fit(ZScoreTransform, Y; dims=1)
    return OutputTransform(transform, output)
end

function (transform::OutputTransform)(df::DataFrame)
    Y = _dataframe_to_array(df, transform.output)
    return StatsBase.transform(transform.transform, Y)
end

function inverse_transform(Y::Array, transform::OutputTransform)
    return StatsBase.reconstruct(transform.transform, Y)
end