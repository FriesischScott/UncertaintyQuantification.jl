struct FORM
    n::Integer
    tol::Real

    function FORM(n::Integer=10, tol::Real=1e-3)
        return new(n, tol)
    end
end

function probability_of_failure(
    models::Union{Array{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Array{<:UQInput},UQInput},
    sim::FORM,
)
    # create reference point in standard normal space origin
    random_names = names(filter(i -> isa(i, RandomUQInput), inputs))
    y::Vector{Float64} = zeros(length(random_names))

    G = [models..., Model(performance, :performance)]

    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)
    parameters =
        !isempty(deterministic_inputs) ? sample(deterministic_inputs, 1) : DataFrame()

    β::Float64 = 0.0
    h₀::Float64 = Inf

    # compute pf through HLRF algorithm
    for it in 1:(sim.n)
        physical = DataFrame(random_names .=> y)
        to_physical_space!(inputs, physical)
        physical = hcat(physical, parameters)

        H = gradient_in_standard_normal_space(G, inputs, physical, :performance)

        H = map(n -> H[n], random_names)

        evaluate!(G, physical)

        h = physical[1, :performance]

        α = H ./ norm(H)

        if it == 1
            h₀ = h
        end

        β = norm(y)

        y_next = -α * (β + (h / norm(H)))
        Δ = norm(y - y_next)

        if Δ < sim.tol
            break
        end

        y = y_next
    end

    pf = cdf(Normal(), (h₀ > 0 ? -1 : 1) * β)
    dp = DataFrame(random_names .=> y)
    to_physical_space!(inputs, dp)

    return pf, β, dp
end
