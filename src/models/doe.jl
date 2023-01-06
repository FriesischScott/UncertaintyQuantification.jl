struct BoxBehnken <: AbstractDesignOfExperiments
end

struct CentralCompositeFaceCentered <: AbstractDesignOfExperiments
end

#inscribed (values of the underlying two lvl factorial will be 1/sqrt(2), not 1)
struct CentralComposite <: AbstractDesignOfExperiments
end

struct FullFactorial <: AbstractDesignOfExperiments
    levels::Array{Integer}
    function FullFactorial(levels)
        return if length(levels) > 0
            new(levels)
        else
            error("levels must hold at least one element")
        end
    end
end

struct TwoLevelFactorial <: AbstractDesignOfExperiments
end

function full_factorial_matrix!(
    m::Matrix, levels::Array, pt_index::Int64, var_index::Int64, block::Int64
)
    if (var_index > length(levels))
        return nothing
    end

    #block size whith same var value
    new_block = (Int64)(block / levels[var_index])

    for i in 1:levels[var_index]
        for j in 1:new_block
            m[(Int64)(new_block * (i - 1) + j + pt_index), var_index] =
                (i - 1) / (levels[var_index] - 1)
        end
        full_factorial_matrix!(
            m, levels, (Int64)(pt_index + new_block * (i - 1)), var_index + 1, new_block
        )
    end
end

function binary_matrix(n::Int64)
    m = ones(2^n, n)

    for i in 1:(2^n)
        for j in 1:n
            m[i, j] = (i - 1) % (2^(n - j + 1)) < 2^(n - j)
        end
    end

    return m
end

function sample(inputs::Array{<:UQInput}, design::AbstractDesignOfExperiments)
    random_inputs = filter(i -> isa(i, RandomUQInput), inputs)

    deterministic_inputs = filter(i -> isa(i, DeterministicUQInput), inputs)

    n_rv = mapreduce(dimensions, +, random_inputs)

    u = doe_samples(design, n_rv)

    samples = quantile.(Normal(), u)
    samples = DataFrame(names(random_inputs) .=> eachcol(samples))

    if !isempty(deterministic_inputs)
        #add deterministic inputs to each row (each point)
        samples = hcat(samples, sample(deterministic_inputs, size(samples, 1)))
    end

    to_physical_space!(inputs, samples)

    return samples
end

function doe_samples(design::FullFactorial, inputs_count::Int64)
    pt_count = reduce(*, design.levels)

    points = ones(pt_count, length(design.levels))

    full_factorial_matrix!(points, design.levels, 0, 1, pt_count)

    return points
end

function doe_samples(design::TwoLevelFactorial, inputs_count::Int64)
    return binary_matrix(inputs_count)
end

function doe_samples(design::CentralCompositeFaceCentered, inputs_count::Int64)
    two_lvl_points = binary_matrix(inputs_count)
    axial_points = ones(2 * inputs_count, inputs_count) .* 0.5

    for i in 1:inputs_count
        for j in 1:2
            axial_points[2 * i + j - 2, i] = j % 2
        end
    end

    return vcat(two_lvl_points, axial_points)
end

function doe_samples(design::CentralComposite, inputs_count::Int64)
    two_lvl_points = binary_matrix(inputs_count)
    axial_points = ones(2 * inputs_count, inputs_count) .* 0.5

    for i in eachrow(two_lvl_points)
        for j in eachindex(i)
            i[j] = 0.5 + 1 / (2.0 * sqrt(2)) * (-1)^(1 + (i[j]))
        end
    end

    for i in 1:inputs_count
        for j in 1:2
            axial_points[2 * i + j - 2, i] = j % 2
        end
    end

    return vcat(two_lvl_points, axial_points)
end

function doe_samples(design::BoxBehnken, inputs_count::Int64)
end
