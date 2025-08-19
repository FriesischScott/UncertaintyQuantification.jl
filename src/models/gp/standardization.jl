"""
    AbstractInputTransform

Abstract type for input transformations used to describe how input features (columns of a DataFrame) should be
preprocessed before fitting a model (e.g. no transform, z-score standardization).
"""
abstract type AbstractInputTransform end

"""
    AbstractOutputTransform

Abstract type for output transformations used to describe how model output (columns of a DataFrame) should be
preprocessed before fitting a model (e.g. no transform, z-score standardization).
"""
abstract type AbstractOutputTransform end

"""
    DataStandardization(input::AbstractInputTransform, output::AbstractOutputTransform)

Container that holds the input and output transformation strategies to be applied
to a dataset.
"""
struct DataStandardization
    input_transform::AbstractInputTransform
    output_transform::AbstractOutputTransform
end

DataStandardization() = DataStandardization(NoInputTransform(), NoOutputTransform())

# ---------------- Input transforms ----------------

"""
    NoInputTransform(input)

No transformation is applied to the specified input columns.
"""
struct NoInputTransform <: AbstractInputTransform end

"""
    ZScoreInputTransform(input)

Applies z-score standardization (mean 0, variance 1) to the specified input columns.
"""

struct ZScoreInputTransform <: AbstractInputTransform end

# ---------------- Output transforms ----------------

"""
    NoOutputTransform(output)

No transformation is applied to the specified output column.
"""
struct NoOutputTransform <: AbstractOutputTransform end

"""
    ZScoreOutputTransform(output)

Applies z-score standardization (mean 0, variance 1) to the specified output column.
Provides both forward (`f`) and inverse (`f⁻¹`) transformations.
"""
struct ZScoreOutputTransform <: AbstractOutputTransform end

# ---------------- Builders ----------------

"""
    build_datatransform(data::DataFrame, transform::AbstractInputTransform)

Builds a transformation function `f(df::DataFrame) -> Array` that converts
the specified input columns of `df` into an array, optionally applying a
standardization transform.
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

"""
    build_datatransform(data::DataFrame, transform::AbstractOutputTransform)

Builds a tuple `(f, f⁻¹)` of transformation functions for the specified output column:

- `f(df)` applies the output transformation to data.
- `f⁻¹(Y)` reverses the transformation for predictions.
"""
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