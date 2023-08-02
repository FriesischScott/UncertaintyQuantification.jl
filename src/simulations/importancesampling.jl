struct ImportanceSampling
    sampling::AbstractMonteCarlo
end

# function sample(inputs::Vector{<:UQInput}, models::Union{Vector{<:UQModel},UQModel}, performance::Function, design_point_from::FORM, sim::ImportanceSampling)
#     # compute design point in standard normal space dp and important direction α
#     pf, β, dp, α = probability_of_failure(models, performance, inputs, design_point_from)
#     dp = DataFrame([dp])
#     to_standard_normal_space!(inputs, dp)
#     α *= -1

#     # generate [n x m] samples in SNS
#     random_names = names(filter(i -> isa(i, RandomUQInput), inputs))
#     sns_inputs = [RandomVariable(Normal(), random_name) for random_name in random_names]
#     Z = Matrix(sample(sns_inputs, sim.sampling))

#     # samples perpendicular to important direction
#     Z -= (Z * α) * α'

#     # force in direction parallel to important direction
#     b = exp(-(β^2/2)) / (cdf(Normal(), -β) * sqrt(2 * π))
#     v = 2(b - β) # no clue how to choose c, here c = 2
#     Zf = Matrix(sample([RandomVariable(Normal(b, v), :Zf)], sim.sampling))

#     Z += Zf * α'

#     weights = DataFrame(ones(sim.sampling.n, 3), [:f, :h, :w]) 
#     weights.f .= pdf.(Normal(), Zf)
#     weights.h .= pdf.(Normal(), (Zf .- b)/v)
#     weights.w = weights.f ./ weights.h

#     random_samples = DataFrame(Z, random_names)
#     to_physical_space!(filter(i -> isa(i, RandomUQInput), inputs), random_samples)
#     deterministic_samples = sample(filter(i -> isa(i, DeterministicUQInput), inputs), sim.sampling)

#     samples = hcat(deterministic_samples, random_samples)

#     return samples, weights
# end

function sample(inputs::Vector{<:UQInput}, proposals::Vector{<:UQInput}, sim::ImportanceSampling)
    # sample from proposal distributions
    samples = sample(proposals, sim.sampling)

    # calculate weights
    weights = DataFrame(ones(size(samples, 1)..., 3), [:f, :h, :w]) 
    weights.f = prod.(eachrow(pdf(inputs, samples)))
    weights.h = prod.(eachrow(pdf(proposals, samples)))
    weights.w = weights.f ./ weights.h

    return samples, weights
end

function proposaldistribution(inputs::Vector{<:UQInput}, models::Union{Vector{<:UQModel},UQModel}, performance::Function, sim::FORM)
    pf, β, dp = probability_of_failure(models, performance, inputs, sim)

    proposals = empty(inputs)

    for input in inputs

        if isa(input, RandomVariable)
            push!(proposals, RandomVariable(Normal(dp[input.name], √(var(input))), input.name))
            continue
        end
        if isa(input, JointDistribution)
            # not sure how to handle this, could create copula of same type with gaussian marginals or independent gaussian distributions
            # create independent gaussian distributions from marginals...
            [push!(proposals, RandomVariable(Normal(dp[marginal.name], √(var(marginal))), marginal.name)) for marginal in input.marginals]
            continue
        end

        push!(proposals, input)
    end

    return proposals
end

proposaldistribution(inputs::Vector{<:UQInput}, performance::Function, sim::FORM) = proposaldistribution(inputs, UQModel[], performance, sim)
proposaldistribution(input::UQInput, models::Union{Vector{<:UQModel},UQModel}, performance::Function, sim::FORM) = proposaldistribution([input], models, performance, sim)
proposaldistribution(input::UQInput, performance::Function, sim::FORM) = proposaldistribution([input], UQModel[], performance, sim) 

# function proposaldistribution(inputs::Vector{<:UQInput}, user_proposals::Vector{<:UQInput})
#     user_proposals_names = names(user_proposals)
#     proposals = empty(inputs)

#     for input in inputs

#         if input.name in user_proposals_names
#             i = findfirst(isequal(input.name), user_proposals_names)
#             push!(proposals, user_proposals[i])
#             continue
#         end

#         push!(proposals, input)
#     end

#     return proposals
# end

# proposaldistribution(input::UQInput, user_proposals::Vector{<:UQInput}) = proposaldistribution([input], user_proposals) 
# proposaldistribution(inputs::Vector{<:UQInput}, user_proposal::UQInput) = proposaldistribution(inputs, [user_proposal]) 
# proposaldistribution(input::UQInput, user_proposal::UQInput) = proposaldistribution([input], [user_proposal]) 
