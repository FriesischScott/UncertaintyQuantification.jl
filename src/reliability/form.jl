struct FORM <: AbstractSimulation
    n::Integer
    tol::Real
    fdm::FiniteDifferencesMethod

    function FORM(
        n::Integer=10,
        tol::Real=1e-3;
        fdm::FiniteDifferencesMethod=CentralFiniteDifferences(3),
    )
        return new(n, tol, fdm)
    end
end

function probability_of_failure(
    models::Union{Vector{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Vector{<:UQInput},UQInput},
    sim::FORM,
)
    models, inputs = wrap.([models, inputs])

    # create reference point in standard normal space origin
    random_names = names(filter(i -> (isa(i, RandomUQInput) || isa(i, ProbabilityBox)), inputs))
    y::Vector{Float64} = zeros(length(random_names))

    deterministic_inputs = filter(i -> (isa(i, DeterministicUQInput) || isa(i, Interval)), inputs)
    parameters =
        !isempty(deterministic_inputs) ? sample(deterministic_inputs, 1) : DataFrame()

    α = Vector{Float64}(undef, length(random_names))
    β::Float64 = 0.0
    h₀::Float64 = Inf

    # compute pf through HLRF algorithm
    for it in 1:(sim.n)
        physical = DataFrame(random_names .=> y)
        physical = hcat(physical, parameters)

        H = gradient_in_standard_normal_space(
            models, inputs, physical, performance, fdm=sim.fdm
        )

        H = map(n -> H[n], random_names)

        to_physical_space!(inputs, physical)

        evaluate!(models, physical)
        h = performance(physical)[1]

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

    dp = (; zip(random_names, collect(dp[1, :]))...)
    α = (; zip(random_names, collect(α))...)
    return pf, β, dp, α
end
