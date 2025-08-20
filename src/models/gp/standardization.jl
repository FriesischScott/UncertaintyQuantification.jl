"""
Input/output transformations for datasets.

- `AbstractInputTransform` / `AbstractOutputTransform`: base types for input and output preprocessing.
- `DataStandardization`: holds the chosen input and output transformations.
- `build_datatransform(data, input/output, transform)`: returns functions that apply (and, for outputs, invert) the transformations to a `DataFrame`.

Predefined transforms include:
    - `NoInputTransform` / `NoOutputTransform`: no change.
    - `ZScoreInputTransform` / `ZScoreOutputTransform`: standardize to zero mean, unit variance.
"""

abstract type AbstractInputTransform end
abstract type AbstractOutputTransform end

struct DataStandardization
    input_transform::AbstractInputTransform
    output_transform::AbstractOutputTransform
end

DataStandardization() = DataStandardization(NoInputTransform(), NoOutputTransform())

# ---------------- Input transforms ----------------

struct NoInputTransform <: AbstractInputTransform end
struct ZScoreInputTransform <: AbstractInputTransform end

# ---------------- Output transforms ----------------

struct NoOutputTransform <: AbstractOutputTransform end
struct ZScoreOutputTransform <: AbstractOutputTransform end

# ---------------- Builders ----------------

"""
build_datatransform(data, input/output, transform)

Returns a function (or pair of functions for outputs) that applies the specified transformation to the dataset.
"""
function build_datatransform(
    data::DataFrame, 
    input::Union{Symbol, Vector{<:Symbol}}, 
    transform::NoInputTransform
)
    f(df::DataFrame) = _dataframe_to_array(df, input)
    return f
end

build_datatransform(
    data::DataFrame, 
    input::Union{UQInput, Vector{<:UQInput}}, 
    transform::NoInputTransform
 ) = build_datatransform(data, names(input), transform)

function build_datatransform(
    data::DataFrame,
    input::Union{Symbol, Vector{<:Symbol}},
    transform::ZScoreInputTransform
)
    input_array = _dataframe_to_array(data, input)
    zscore_transform = fit(ZScoreTransform, input_array; dims=1)
    f(df::DataFrame) = StatsBase.transform(
            zscore_transform, 
            _dataframe_to_array(df, input)
        )
    return f
end

build_datatransform(
    data::DataFrame, 
    input::Union{UQInput, Vector{<:UQInput}}, 
    transform::ZScoreInputTransform
 ) = build_datatransform(data, names(input), transform)

function build_datatransform(
    data::DataFrame, 
    output::Symbol, 
    transform::NoOutputTransform
)
    f(df::DataFrame) = _dataframe_to_array(df, output)
    f⁻¹(Y::AbstractArray) = Y
    return (f, f⁻¹)
end

function build_datatransform(
    data::DataFrame, 
    output::Symbol,
    transform::ZScoreOutputTransform
)
    output_array = _dataframe_to_array(data, output) # will fail if Y is not an array
    zscore_transform = fit(ZScoreTransform, output_array; dims=1)
    f(df::DataFrame) = StatsBase.transform(
            zscore_transform, 
            _dataframe_to_array(df, output)
        )
    f⁻¹(Y::AbstractArray) = StatsBase.reconstruct(zscore_transform, Y)
    return return (f, f⁻¹)
end

function build_datatransforms(
    data::DataFrame,
    input::Union{Symbol, Vector{<:Symbol}, UQInput, Vector{<:UQInput}},
    output::Symbol,
    ds::DataStandardization
)
    fᵢ = build_datatransform(data, input, ds.input_transform)
    fₒ, fₒ⁻¹ = build_datatransform(data, output, ds.output_transform)
    return (fᵢ, fₒ, fₒ⁻¹)
end


# ### Standard normal transform
# struct UQInputTransform <: AbstractInputTransform
#     uqinput::Union{UQInput, Vector{<:UQInput}}
# end

# function (transform::UQInputTransform)(df::DataFrame)
#     df_copy = copy(df)
#     uqinput_names = names(transform.uqinput)
#     to_standard_normal_space!(transform.uqinput, df_copy)
#     # X is a Matrix for multiple inputs, else it is a Vector 
#     X = _dataframe_to_array(df_copy, uqinput_names)
#     return X
# end

# ---------------- Utility ----------------

_dataframe_to_array(df::DataFrame, name::Symbol) = df[:, name]

function _dataframe_to_array(df::DataFrame, names::Vector{<:Symbol})
    length(names) == 1 ? x = _dataframe_to_array(df, only(names)) : x = RowVecs(Matrix(df[:, names]))
    return x
end