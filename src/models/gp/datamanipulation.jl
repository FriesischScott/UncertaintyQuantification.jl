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
    output::Symbol,
    normalize_inp::Bool=false,
    normalize_out::Bool=false
)
    return _handle_gp_input(
        data, inputs, [output],
        normalize_inp, normalize_out
        )
end