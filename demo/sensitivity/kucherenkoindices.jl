include("../../src/UncertaintyQuantification.jl")
using .UncertaintyQuantification
using LinearAlgebra

# Testing Kucherenko indices - Test Case 1: Linear model with correlated variables from Kucherenko et al. (2012)

ρ = 0.5  # correlation coefficient
σ = 2.0  # standard deviation for x3


Σ = [1.0  0.0  0.0;
     0.0  1.0  ρ*σ;
     0.0  ρ*σ   σ^2]

D = sqrt.(diag(Σ)) 
R = Σ ./ (D * D') 

marginals = UncertaintyQuantification.RandomVariable[
    UncertaintyQuantification.RandomVariable(Normal(0, 1), :x1),    # std = 1
    UncertaintyQuantification.RandomVariable(Normal(0, 1), :x2),    # std = 1  
    UncertaintyQuantification.RandomVariable(Normal(0, σ), :x3)     # std = σ
]

inputs = [
    UncertaintyQuantification.JointDistribution(marginals, UncertaintyQuantification.GaussianCopula(R))
]

model = UncertaintyQuantification.Model(df -> df.x1 .+ df.x2 .+ df.x3, :y)

sim = UncertaintyQuantification.MonteCarlo(500000)

try
    indices = UncertaintyQuantification.kucherenkoindices([model], inputs, [:y], sim, num_bins=60)

    println(indices)

    # Analytical calculations
    denom = 2 + σ^2 + 2*ρ*σ
    S1_analytical = 1 / denom
    S2_analytical = (1 + ρ*σ)^2 / denom
    S3_analytical = (σ + ρ)^2 / denom
    ST1_analytical = 1 / denom
    ST2_analytical = (1 - ρ^2) / denom
    ST3_analytical = σ^2 * (1 - ρ^2) / denom
    

    computed_values = Dict()
    for i in 1:size(indices, 1)
        var_name = string(indices[i, :Variables])
        computed_values[var_name] = (indices[i, :FirstOrder], indices[i, :TotalEffect])
    end
    
    tol = 0.01
    @assert abs(computed_values["x1"][1] - S1_analytical) < tol "S1 error too large"
    @assert abs(computed_values["x2"][1] - S2_analytical) < tol "S2 error too large"
    @assert abs(computed_values["x3"][1] - S3_analytical) < tol "S3 error too large"
    @assert abs(computed_values["x1"][2] - ST1_analytical) < tol "ST1 error too large"
    @assert abs(computed_values["x2"][2] - ST2_analytical) < tol "ST2 error too large"
    @assert abs(computed_values["x3"][2] - ST3_analytical) < tol "ST3 error too large"
    
    println("Kucherenko indices validation passed - all values within tolerance.")
    
catch e
    println("Error computing Kucherenko indices: $e")
    rethrow(e)
end


