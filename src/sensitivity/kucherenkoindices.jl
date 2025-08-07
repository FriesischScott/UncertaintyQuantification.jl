const _kucherenko_table_types = [Symbol[], Float64[], Float64[]]
const _kucherenko_table_header = [
    "Variables", "FirstOrder", "TotalEffect"
]

"""
- `X::Matrix`: Input sample matrix where each row is a sample point and each column is a variable
- `Y::Vector`: Output sample vector corresponding to the input samples
- `num_bins::Int=10`: Number of bins to use for conditioning
"""
function kucherenkoindices(X::Matrix, Y::Vector, num_bins::Int=10)
    n_samples, n_vars = size(X)
    
    if length(Y) != n_samples
        throw(ArgumentError("Number of output samples must match number of input samples"))
    end
    
    indices = DataFrame(_kucherenko_table_types, _kucherenko_table_header)
    total_var = var(Y)
    
    for i in 1:n_vars
        S_i = _compute_first_order_kucherenko(X, Y, i, num_bins, total_var)
        ST_i = _compute_total_effect_kucherenko(X, Y, i, num_bins, total_var)
        
        push!(indices, [Symbol("X$i"), S_i, ST_i])
    end
    
    return indices
end

function _compute_first_order_kucherenko(X::Matrix, Y::Vector, var_idx::Int, num_bins::Int, total_var::Float64)
    n_samples = length(Y)
    
    x_var = X[:, var_idx]
    
    x_min, x_max = extrema(x_var)
    if x_min ≈ x_max
        return 0.0 
    end
    
    bin_edges = range(x_min, x_max, length=num_bins+1)
    
    bin_means = Float64[]
    bin_weights = Float64[]
    
    for j in 1:num_bins
        if j == num_bins
            mask = (x_var .>= bin_edges[j]) .& (x_var .<= bin_edges[j+1])
        else
            mask = (x_var .>= bin_edges[j]) .& (x_var .< bin_edges[j+1])
        end
        
        n_bin = sum(mask)
        
        if n_bin > 0
            bin_mean = mean(Y[mask])
            push!(bin_means, bin_mean)
            push!(bin_weights, n_bin / n_samples)
        end
    end
    
    if isempty(bin_means)
        return 0.0
    end
    
    weighted_mean = sum(bin_weights .* bin_means)
    
    S_i = sum(bin_weights .* (bin_means .- weighted_mean).^2) / total_var
    
    return S_i
end

function _compute_total_effect_kucherenko(X::Matrix, Y::Vector, var_idx::Int, num_bins::Int, total_var::Float64)
    n_samples, n_vars = size(X)
    
    if n_vars == 1
        return _compute_first_order_kucherenko(X, Y, var_idx, num_bins, total_var)
    end
    
    other_vars = setdiff(1:n_vars, var_idx)
    X_other = X[:, other_vars]
    
    bin_assignments = _assign_multidimensional_bins(X_other, num_bins)
    
    weighted_var_sum = 0.0
    total_weight = 0.0
    
    for bin_id in unique(bin_assignments)
        mask = bin_assignments .== bin_id
        n_bin = sum(mask)
        
        if n_bin > 1
            bin_var = var(Y[mask])
            weight = n_bin / n_samples
            weighted_var_sum += weight * bin_var
            total_weight += weight
        end
    end
    
    ST_i = total_weight > 0 ? weighted_var_sum / total_var : 0.0
    
    return ST_i
end

function _assign_multidimensional_bins(X::Matrix, num_bins::Int)
    n_samples, n_dims = size(X)
    
    X_norm = similar(X)
    for j in 1:n_dims
        x_min, x_max = extrema(X[:, j])
        if x_min ≈ x_max
            X_norm[:, j] .= 0.5
        else
            X_norm[:, j] = (X[:, j] .- x_min) ./ (x_max - x_min)
        end
    end
    
    bin_assignments = zeros(Int, n_samples)
    
    for i in 1:n_samples
        bin_id = 0
        for j in 1:n_dims
            bin_idx = min(floor(Int, X_norm[i, j] * num_bins), num_bins - 1)
            bin_id += bin_idx * (num_bins ^ (j - 1))
        end
        bin_assignments[i] = bin_id + 1
    end
    
    return bin_assignments
end

"""
- `models::Vector{<:UQModel}`: Vector of UQ models to evaluate
- `inputs::Vector{<:UQInput}`: Vector of input distributions
- `outputs::Vector{Symbol}`: Vector of output quantity names
- `sim::AbstractMonteCarlo`: Monte Carlo simulation parameters
- `num_bins::Int=10`: Number of bins for conditioning (default: 10)
"""
function kucherenkoindices(
    models::Vector{<:UQModel},
    inputs::Vector{<:UQInput},
    outputs::Vector{Symbol},
    sim::AbstractMonteCarlo;
    num_bins::Int=10
)
    samples = sample(inputs, sim)
    
    evaluate!(models, samples)
    
    random_names = names(filter(i -> isa(i, RandomUQInput), inputs))
    
    indices = Dict([
        (name, DataFrame(_kucherenko_table_types, _kucherenko_table_header)) for name in outputs
    ])
    
    X = Matrix(samples[:, random_names])
    
    for output in outputs
        Y = Vector(samples[:, output])
        
        output_indices = kucherenkoindices(X, Y, num_bins)
        
        for (i, var_name) in enumerate(random_names)
            output_indices[i, :Variables] = var_name
        end
        
        indices[output] = output_indices
    end
    
    return length(outputs) > 1 ? indices : indices[outputs[1]]
end

function kucherenkoindices(
    models::Vector{<:UQModel},
    inputs::UQInput,
    outputs::Vector{Symbol},
    sim::AbstractMonteCarlo;
    num_bins::Int=10
)
    return kucherenkoindices(models, [inputs], outputs, sim; num_bins=num_bins)
end

function kucherenkoindices(
    models::Vector{<:UQModel},
    inputs::Vector{<:UQInput},
    outputs::Symbol,
    sim::AbstractMonteCarlo;
    num_bins::Int=10
)
    return kucherenkoindices(models, inputs, [outputs], sim; num_bins=num_bins)
end

function kucherenkoindices(
    models::UQModel,
    inputs::Vector{<:UQInput},
    outputs::Symbol,
    sim::AbstractMonteCarlo;
    num_bins::Int=10
)
    return kucherenkoindices([models], inputs, [outputs], sim; num_bins=num_bins)
end

function kucherenkoindices(
    models::Vector{<:UQModel},
    inputs::UQInput,
    outputs::Symbol,
    sim::AbstractMonteCarlo;
    num_bins::Int=10
)
    return kucherenkoindices(models, [inputs], [outputs], sim; num_bins=num_bins)
end

function kucherenkoindices(
    models::UQModel,
    inputs::UQInput,
    outputs::Symbol,
    sim::AbstractMonteCarlo;
    num_bins::Int=10
)
    return kucherenkoindices([models], [inputs], [outputs], sim; num_bins=num_bins)
end
