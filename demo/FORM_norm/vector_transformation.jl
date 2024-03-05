using LinearAlgebra

function gradients_l1_from_l2(gradients_l2)
    # Step 1: Normalize the vector in L2 norm
    norm_l2 = norm(gradients_l2)
    gradients_normalized = gradients_l2 / norm_l2
    
    # Step 2: Scale the vector to L1 norm
    norm_l1 = sum(abs.(gradients_l2))
    gradients_l1 = norm_l1 * gradients_normalized
    
    return gradients_l1
end

# Example usage
gradients_l2 = [1, -2, 3]  # Example gradients in L2 norm
gradients_l1 = gradients_l1_from_l2(gradients_l2)
println("Gradients in L1 norm: ", gradients_l1)


function transform_l2_to_l1(vector_l2, l1_norm)
    # Step 1: Sort the absolute values of the vector components
    sorted_abs_values = sort(abs.(vector_l2), rev=true)
    
    # Step 2: Reconstruct the vector with sorted absolute values
    sorted_vector = sign.(vector_l2) .* sorted_abs_values
    
    # Step 3: Normalize the vector to the desired L1 norm
    norm_l2 = norm(vector_l2)
    norm_ratio = l1_norm / norm_l2
    vector_l1 = norm_ratio * sorted_vector
    
    return vector_l1
end

# Example usage
vector_l2 = [1, -2]  # Example vector in L2 norm
l1_norm = 1  # Desired L1 norm
vector_l1 = transform_l2_to_l1(vector_l2, l1_norm)
println("Vector in L1 norm: ", vector_l1)

using Plots

N_samples = 1000
l2_samples = randn(N_samples, 2)

l1_samples = [transform_l2_to_l1(l2_samples[i,:], l1_norm) for i in 1:N_samples]
l1_samples = transpose(reduce(hcat, l1_samples))

scatter(l2_samples[:,1], l2_samples[:, 2])
scatter!(l1_samples[:,1], l1_samples[:, 2])