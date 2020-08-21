function probability_of_failure(
    models::Union{Array{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Array{<:UQInput},UQInput},
    sim::MonteCarlo,
)

    samples = sample(inputs, sim.n)

    # Models
    for m in models
        evaluate!(m, samples)
    end

    # Probability of failure
    pf = sum(performance(samples) .< 0) / sim.n

    return pf, samples
end

function probability_of_failure(
    models::Union{Array{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Array{<:UQInput},UQInput},
    sim::LineSampling,
)

    if isempty(sim.direction)
        sim.direction = gradient_in_standard_normal_space(
        [models..., Model(x -> -1 * performance(x), :performance)],
        inputs,
        mean(inputs),
        :performance)
    end

    samples = sample(inputs, sim)

    for m in models
        evaluate!(m, samples)
    end

    p = reshape(performance(samples), length(sim.points), sim.lines)

    ϕ = Normal()
    pf = 0
    roots_found = 0
    x = median(sim.points)
    for i = 1:sim.lines
        if all(p[:, i] .< 0)
            @warn "All samples for line $i are inside the failure domain"
            continue
        elseif all(p[:, i] .> 0)
            @warn "All samples for line $i are outside the failure domain"
            continue
        end
        spl = Spline1D(sim.points, p[:, i])
        try
            root = Dierckx.roots(spl)[1]
            pf += cdf.(ϕ, -1 * root)
            roots_found += 1
        catch e
            @warn "Intersection with failure domain not found for line $i ($e)"
        end
    end

    pf /= roots_found

    return pf, samples
end

const Φ = Normal()

function FORM(models::Union{Array{<:UQModel},UQModel},
    performance::Function,
    inputs::Union{Array{<:UQInput},UQInput})

    # create reference point in standard normal space origin
    random_names = filter(i -> isa(i, RandomUQInput), inputs) |> names
    u = DataFrame(zeros(1, length(random_names)))
    rename!(u, random_names)

    physical = to_physical_space(inputs, u)
    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)

    if !isempty(deterministic_inputs)
        physical = hcat(physical, sample(deterministic_inputs, 1))
    end
    
    iteration = 0
    delta = Inf
    beta_last = Inf
    performance_at_origin = Inf

    # compute pf through HLRF algorithm
    while delta > 1e-2 && iteration < 10
        iteration += 1

        # @show u
        @show physical

        gradient = gradient_in_standard_normal_space(
            [models..., Model(performance, :performance)],
            inputs,
            physical,
            :performance) |> Matrix

        for m in models
            evaluate!(m, physical)
        end

        α = gradient / norm(gradient)
        @show α

        b = (Matrix(u) * transpose(α)) * α
        @show b
        u[:,:] = b - physical[1, last(models).name] / norm(gradient) * α

        physical[:, names(u)] = u;

        if iteration == 1
            beta_last = norm(Matrix(u))
            performance_at_origin = physical[1, last(models).name]
        else
            beta = norm(Matrix(u))
            delta = abs((beta - beta_last) / beta_last)
            beta_last = beta
        end

    end

    @show norm(Matrix(u))

    pf = cdf(Φ , (performance_at_origin > 0 ? -1 : 1) * norm(Matrix(u)))
    designPoint = to_physical_space(inputs, u)
    
    return pf, designPoint

end
