"""
Input/output transformations for datasets.

- `AbstractInputTransform` / `AbstractOutputTransform`: base types for input and output preprocessing.
- `DataTransforms`: holds the chosen input and output transformations.
- `build_datatransform(data, input/output, transform)`: returns functions that apply (and, for outputs, invert) the transformations to a `DataFrame`.

Predefined transforms include:
    - `NoInputTransform` / `NoOutputTransform`: no change.
    - `ZScoreInputTransform` / `ZScoreOutputTransform`: standardize to zero mean, unit variance.
"""
abstract type AbstractDataTransform end

# ---------------- Input/Output transforms ----------------
struct IdentityTransform <: AbstractDataTransform end
struct ZScoreTransform <: AbstractDataTransform end
struct UnitRangeTransform <: AbstractDataTransform end
struct StandardNormalTransform <: AbstractDataTransform end

struct InputTransform{T <: AbstractDataTransform} end
InputTransform(::Type{T}) where {T <: AbstractDataTransform} = InputTransform{T}()
InputTransform(x::AbstractDataTransform) = InputTransform(typeof(x))

struct OutputTransform{T <: AbstractDataTransform} end
OutputTransform(::Type{T}) where {T <: AbstractDataTransform} = OutputTransform{T}()
OutputTransform(x::AbstractDataTransform) = OutputTransform(typeof(x))

# ---------------- Struct for bundled transform functions ----------------
struct DataStandardizer
    fᵢ::Function
    fₒ::Function
    fₒ⁻¹::Function
    var_fₒ⁻¹::Function
end

# ---------------- Constructor ----------------
function DataStandardizer(
    data::DataFrame,
    input::Union{Symbol, Vector{<:Symbol}, UQInput, Vector{<:UQInput}},
    output::Symbol,
    input_transform::InputTransform, 
    output_transform::OutputTransform
)
    fᵢ = build_datatransform(data, input, input_transform)
    fₒ, fₒ⁻¹, var_fₒ⁻¹ = build_datatransform(data, output, output_transform)
    return DataStandardizer(fᵢ, fₒ, fₒ⁻¹, var_fₒ⁻¹)
end

# ---------------- Transform builders ----------------
"""
build_datatransform(data, input/output, transform)

Returns a function (or pair of functions for outputs) that applies the specified transformation to a dataframe.
"""
# ---------------- Input ----------------
# No input transformation
function build_datatransform(
    ::DataFrame, 
    input::Union{Symbol, Vector{<:Symbol}}, 
    ::InputTransform{IdentityTransform}
)
    f(df::DataFrame) = to_gp_format(
        dataframe_to_array(df, input)
    )
    return f
end

build_datatransform(
    data::DataFrame, 
    input::Union{UQInput, Vector{<:UQInput}}, 
    transform::InputTransform{IdentityTransform}
 ) = build_datatransform(data, names(input), transform)

 # ZScore input transformation
function build_datatransform(
    data::DataFrame,
    input::Union{Symbol, Vector{<:Symbol}},
    ::InputTransform{ZScoreTransform}
)
    zscore_transform = fit(
        StatsBase.ZScoreTransform, 
        dataframe_to_array(data, input); 
        dims=1
    )
    f(df::DataFrame) = to_gp_format(
        StatsBase.transform(
            zscore_transform, 
            dataframe_to_array(df, input)
        )
    )
    return f
end

build_datatransform(
    data::DataFrame, 
    input::Union{UQInput, Vector{<:UQInput}}, 
    transform::InputTransform{ZScoreTransform}
 ) = build_datatransform(data, names(input), transform)

# UnitRange input transformation
function build_datatransform(
    data::DataFrame,
    input::Union{Symbol, Vector{<:Symbol}},
    ::InputTransform{UnitRangeTransform}
)
    unitrange_transform = fit(
        StatsBase.UnitRangeTransform, 
        dataframe_to_array(data, input); 
        dims=1
    )
    f(df::DataFrame) = to_gp_format(
        StatsBase.transform(
            unitrange_transform, 
            dataframe_to_array(df, input)
        )
    )
    return f
end

build_datatransform(
    data::DataFrame, 
    input::Union{UQInput, Vector{<:UQInput}}, 
    transform::InputTransform{UnitRangeTransform}
 ) = build_datatransform(data, names(input), transform)

# SNS input transform
function build_datatransform(
    ::DataFrame,
    input::Union{UQInput, Vector{<:UQInput}},
    ::InputTransform{StandardNormalTransform}
)
    function f(df::DataFrame)
        df_copy = copy(df)
        to_standard_normal_space!(input, df_copy)
        return to_gp_format(
            dataframe_to_array(df_copy, names(input))    
        )
    end
    return f
end

# ---------------- Output ----------------
# No output transformation
function build_datatransform(
    ::DataFrame, 
    output::Symbol, 
    ::OutputTransform{IdentityTransform}
)
    f(df::DataFrame) = to_gp_format(
        dataframe_to_array(df, output)
    )
    f⁻¹(Y::AbstractArray) = Y
    var_f⁻¹(Y::AbstractArray) = Y
    return (f, f⁻¹, var_f⁻¹)
end

# ZScore output transformation
function build_datatransform(
    data::DataFrame, 
    output::Symbol,
    ::OutputTransform{ZScoreTransform}
)
    zscore_transform = fit(
        StatsBase.ZScoreTransform, 
        dataframe_to_array(data, output); 
        dims=1
    )
    f(df::DataFrame) = to_gp_format(
        StatsBase.transform(
            zscore_transform, 
            dataframe_to_array(df, output)
        )
    )
    f⁻¹(Y::AbstractArray) = StatsBase.reconstruct(zscore_transform, Y)
    var_f⁻¹(Y::AbstractArray) = only(zscore_transform.scale)^2 * Y 
    return (f, f⁻¹, var_f⁻¹)
end

# UnitRange output transformation
function build_datatransform(
    data::DataFrame, 
    output::Symbol,
    ::OutputTransform{UnitRangeTransform}
)
    unitrange_transform = fit(
        StatsBase.UnitRangeTransform, 
        dataframe_to_array(data, output); 
        dims=1
    )
    f(df::DataFrame) = to_gp_format(
        StatsBase.transform(
            unitrange_transform, 
            dataframe_to_array(df, output)
        )
    )
    f⁻¹(Y::AbstractArray) = StatsBase.reconstruct(unitrange_transform, Y)
    var_f⁻¹(Y::AbstractArray) = only(unitrange_transform.scale)^2 * Y 
    return (f, f⁻¹, var_f⁻¹)
end

function build_datatransform(
    ::DataFrame, 
    ::Symbol, 
    ::OutputTransform{StandardNormalTransform}
)
    throw(ArgumentError(
        "StandardNormalTransform is only valid for input transforms."
    ))
end

# ---------------- Utility ----------------
to_gp_format(x::Vector) = x
to_gp_format(x::Matrix) = RowVecs(x)

dataframe_to_array(df::DataFrame, name::Symbol) = df[:, name]
dataframe_to_array(df::DataFrame, names::Vector{<:Symbol}) = length(names) == 1 ? x = dataframe_to_array(df, only(names)) : x = Matrix(df[:, names])