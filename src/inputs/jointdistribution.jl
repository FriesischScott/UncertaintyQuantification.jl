struct JointDistribution <: RandomUQInput
    marginals::Vector{RandomVariable}
    copula::Copula

    function JointDistribution(marginals::Vector{RandomVariable}, copula::Copula)
        length(marginals) == dimensions(copula) ||
            error("Dimension mismatch between copula and marginals")

        return new(marginals, copula)
    end
end

function sample(jd::JointDistribution, n::Integer=1)
    u = sample(jd.copula, n)

    samples = DataFrame()

    for (i, rv) in enumerate(jd.marginals)
        samples[!, rv.name] = quantile.(rv.dist, u[:, i])
    end

    return samples
end


function sample_conditional_copula(joint::JointDistribution, var_names::Vector{Symbol}, x_values::Vector{Float64}, N::Int)
    length(var_names) != length(x_values) && error("Number of variable names ($(length(var_names))) must match number of values ($(length(x_values)))")
    !isa(joint.copula, GaussianCopula) && error("This function for conditional sampling is only implemented for Gaussian copulas.")
    
    marginals, copula, R = joint.marginals, joint.copula, joint.copula.correlation
    all_var_names, d = [marginal.name for marginal in marginals], length(marginals)
    
    v_indices = Int[]
    for var_name in var_names
        idx = findfirst(==(var_name), all_var_names)
        idx === nothing && error("Variable $var_name not found in joint distribution. Available variables: $all_var_names")
        push!(v_indices, idx)
    end
    
    length(v_indices) >= d && error("All $(d) variables are fixed - need at least one free variable to sample")
    
    w_indices = setdiff(1:d, v_indices)
    
    z_v = zeros(length(v_indices))
    for (i, (idx, x_val)) in enumerate(zip(v_indices, x_values))
        z_v[i] = quantile(Normal(0,1), cdf(marginals[idx].dist, x_val))
    end
    
    Σ_vv, Σ_wv, Σ_vw, Σ_ww = R[v_indices, v_indices], R[w_indices, v_indices], R[v_indices, w_indices], R[w_indices, w_indices]
    
    if length(v_indices) == 1
        μ_cond = (Σ_wv / Σ_vv[1,1]) * z_v[1]
        Σ_cond = Σ_ww .- (Σ_wv * Σ_vw) / Σ_vv[1,1]
    else
        Σ_vv_inv = inv(Σ_vv)
        μ_cond = Σ_wv * Σ_vv_inv * z_v
        Σ_cond = Σ_ww .- Σ_wv * Σ_vv_inv * Σ_vw
    end
    
    μ_cond = vec(μ_cond)
    Z_w = rand(MvNormal(μ_cond, Symmetric(Σ_cond)), N)'
    
    samples = DataFrame()
    for (i, var_name) in enumerate(var_names)
        samples[!, var_name] = fill(x_values[i], N)
    end
    for (i, idx) in enumerate(w_indices)
        samples[!, all_var_names[idx]] = quantile.(marginals[idx].dist, cdf.(Normal(), Z_w[:,i]))
    end
    
    return samples
end

function sample_conditional_copula(joint::JointDistribution, var_name::Symbol, x_v::Float64, N::Int)
    return sample_conditional_copula(joint, [var_name], [x_v], N)
end

function to_physical_space!(jd::JointDistribution, x::DataFrame)
    correlated_cdf = to_copula_space(jd.copula, Matrix{Float64}(x[:, names(jd)]))
    for (i, rv) in enumerate(jd.marginals)
        x[!, rv.name] = quantile.(rv.dist, correlated_cdf[:, i])
    end
    return nothing
end

function to_standard_normal_space!(jd::JointDistribution, x::DataFrame)
    for rv in jd.marginals
        x[!, rv.name] = cdf.(rv.dist, x[:, rv.name])
    end
    uncorrelated_stdnorm = to_standard_normal_space(
        jd.copula, Matrix{Float64}(x[:, names(jd)])
    )
    for (i, rv) in enumerate(jd.marginals)
        x[!, rv.name] = uncorrelated_stdnorm[:, i]
    end
    return nothing
end

names(jd::JointDistribution) = vec(map(x -> x.name, jd.marginals))

mean(jd::JointDistribution) = mean.(jd.marginals)

dimensions(jd::JointDistribution) = dimensions(jd.copula)
