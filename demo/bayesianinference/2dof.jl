using UncertaintyQuantification

function example2dof(θ::AbstractArray)
    λ1_noise = ((θ[1] + 2 * θ[2]) + sqrt(θ[1]^2 + 4 * θ[2]^2)) / 2 + rand(Normal())
    λ2_noise = ((θ[1] + 2 * θ[2]) - sqrt(θ[1]^2 + 4 * θ[2]^2)) / 2 + rand(Normal(0, 0.5))
    return [λ1_noise, λ2_noise]
end

function likelihood(λ_noisy::AbstractArray, λ_model::AbstractArray)
    for dim in 1:size(λ_model, 2)
        for iobservation in 1:size(λ_noisy, 1)
            #blablablabö
        end
    end

    return [λ1_noise, λ2_noise]
end
