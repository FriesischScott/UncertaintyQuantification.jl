
struct TwoLevelFactorial <: AbstractDesignOfExperiments end

struct FullFactorial <: AbstractDesignOfExperiments
    levels::Vector{<:Integer}
    function FullFactorial(levels::Vector{<:Integer})
        any(levels .< 2) && error("Columns must be >= 2")
        return new(levels)
    end
end

#this is a 2lvl fractional factorial
struct FractionalFactorial <: AbstractDesignOfExperiments
    columns::Vector{String}
end

struct BoxBehnken <: AbstractDesignOfExperiments
    centers::Int
    function BoxBehnken(centers::Int=-1)
        return new(centers)
    end
end

# struct CentralComposite <: AbstractDesignOfExperiments
#     type::String
#     function CentralComposite(type)
#         return if type == "inscribed" || type == "face_centered"
#             new(type)
#         else
#             error("non existing type")
#     end
# end

#@enum CCType face_centered = 1 inscribed = 2

struct CentralComposite <: AbstractDesignOfExperiments
    CCtype::Symbol
    function CentralComposite(CCtype::Symbol)
        CCtype != :inscribed &&
            CCtype != :face &&
            error("CCtype must be :inscribed or :face.")
        return new(CCtype)
    end
end

function full_factorial_matrix(levels::Vector{<:Integer})
    ranges = [range(0.0, 1.0; length=l) for l in levels]
    return mapreduce(t -> [t...]', vcat, Iterators.product(ranges...))
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

function doe_samples(design::FullFactorial, _::Int=0)
    return full_factorial_matrix(design.levels)
end

function doe_samples(_::TwoLevelFactorial, rvs::Int)
    return full_factorial_matrix(ones(Int, rvs) .* 2)
end

function doe_samples(design::FractionalFactorial, _::Int=0)
    vars, varindex, gen = extract_vars_and_generator(design.columns)   #vars holds all variables, gen holds generators(multi-variable strings)

    nvars = length(vars)

    #ff for single vars
    m1 = full_factorial_matrix(ones(Int, nvars) .* 2)
    m2 = zeros(2^nvars, length(gen))
    m2_index = 1

    #add generator columns
    for i in eachindex(gen)
        for row in 1:(2^nvars)
            entry = 1
            for j in eachindex(gen[i])
                letter_index = findfirst(gen[i][j], vars)
                if (m1[row, letter_index] == 0)
                    entry += 1
                end
            end
            m2[row, m2_index] = entry % 2
        end
        m2_index += 1
    end
    #sorting columns to original order
    m = zeros(2^nvars, 1)
    v = 1
    g = 1
    for i in eachindex(design.columns)
        if (findfirst(isequal(i), varindex) === nothing)
            m = hcat(m, m2[:, g])
            g += 1
        else
            m = hcat(m, m1[:, v])
            v += 1
        end
    end
    return m = m[:, 2:end]
end

function extract_vars_and_generator(columns::Vector{String})
    gen = []
    vars = ""
    varindex = Vector{Int}(undef, 0)
    for i in eachindex(columns)
        if (length(columns[i]) == 0)
            error("each String in columns must hold at least one character")
        elseif (length(columns[i]) == 1)
            if (!contains(vars, columns[i][1]))
                vars = vars * columns[i][1]     #collect variables
                push!(varindex, i)
            end
        else
            push!(gen, columns[i])   #collect generators
        end
    end
    for i in eachindex(columns)
        for j in eachindex(columns[i])
            columns[i]
            if !contains(vars, columns[i][j])
                error(
                    "FractionalFactorial: a variable cannot only be used in combination with others, must have its own column.",
                )
            end
        end
    end
    return vars, varindex, gen
end

function doe_samples(design::BoxBehnken, rvs::Int)
    if (design.centers < 0)
        centers = [1 1 3 3 6 6 6 8 10 10 12][min(rvs, 11)]
    else
        centers = design.centers
    end

    ff_columns, block_format = get_block_format(rvs)
    ff = FractionalFactorial(ff_columns)
    ff_m = doe_samples(ff)
    bb = zeros(1, rvs)

    for i in axes(block_format)[1] #iterating the blocks
        m_help = ones(size(ff_m, 1), rvs) * 0.5
        for j in axes(block_format)[2] # iterating over vars of current block
            for k in axes(ff_m)[1] # iterating each row of current block
                m_help[k, block_format[i, j]] = ff_m[k, j]
            end
        end
        bb = vcat(bb, m_help)
    end
    return vcat(bb[2:end, :], ones(centers, rvs) * 0.5)
end

function get_block_format(nvars::Int)
    if nvars == 3
        return ["a", "b"], [1 2; 1 3; 2 3]
    elseif nvars == 4
        return ["a", "b"], [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]
    elseif nvars == 5
        return ["a", "b"], [1 2; 1 3; 1 4; 1 5; 2 3; 2 4; 2 5; 3 4; 3 5; 4 5]
    elseif nvars == 6
        return ["a", "b", "c"], [1 4 5; 1 3 6; 1 2 4; 2 5 6; 2 3 5; 3 4 6]
    elseif nvars == 7
        return ["a", "b", "c"], [1 2 3; 1 4 6; 1 5 7; 2 5 6; 2 4 7; 3 4 5; 3 6 7]
    elseif nvars == 9
        return ["a", "b", "c"],
        [
            1 2 3
            1 4 5
            1 6 7
            1 8 9
            1 4 7
            2 4 6
            2 5 8
            2 7 9
            2 5 8
            3 5 7
            2 4 9
            3 6 8
            3 6 9
            4 7 8
            5 6 9
        ]
    elseif nvars == 10
        return ["a", "b", "c", "d"],
        [
            1 2 5 10
            1 4 7 8
            1 3 6 9
            1 8 9 10
            2 6 7 10
            2 3 7 8
            2 4 6 9
            3 4 5 10
            3 5 7 9
            4 5 6 8
        ]
    elseif nvars == 11
        return ["a", "b", "c", "d", "abcd"],
        [
            1 4 8 9 10
            1 3 6 10 11
            1 2 4 7 11
            1 2 3 5 8
            1 5 6 7 9
            2 5 9 10 11
            2 3 4 6 9
            2 6 7 8 10
            3 7 8 9 11
            3 4 5 7 10
            4 5 6 8 11
        ]
    elseif nvars == 12
        return ["a", "b", "c", "d"],
        [
            1 2 5 7
            1 7 8 11
            1 3 9 10
            1 4 6 12
            2 3 6 8
            2 8 9 12
            2 4 10 11
            3 4 7 9
            3 5 11 12
            4 5 8 10
            5 6 9 11
            6 7 10 12
        ]
    elseif nvars == 16
        return ["a", "b", "c", "d"],
        [
            1 2 6 9
            1 4 5 12
            1 9 10 14
            1 5 8 16
            1 3 13 15
            1 7 11 13
            2 3 7 10
            2 4 14 16
            2 5 6 13
            2 10 11 15
            2 8 12 14
            3 4 8 11
            3 6 7 14
            3 11 12 16
            3 5 9 15
            4 7 8 15
            4 9 12 13
            4 6 10 16
            5 10 13 14
            5 7 9 11
            6 11 14 15
            6 8 10 12
            7 12 15 16
            8 9 13 16
        ]
    elseif nvars > 0
        return ["a", "b"], calc_blocks(nvars)
    end
end

function calc_blocks(nvars::Int)
    out = zeros(Int, 1, 2)
    for i in 1:(nvars - 1)
        for j in (i + 1):nvars
            out = vcat(out, [i j])
        end
    end
    return out[2:end, :]
end

function doe_samples(design::CentralComposite, rvs::Int)
    two_lvl_points = full_factorial_matrix(ones(Int, rvs) .* 2)
    axial_points = ones(2 * rvs + 1, rvs) .* 0.5

    if (design.CCtype == :inscribed)
        for i in eachrow(two_lvl_points)
            for j in eachindex(i)
                #shortening "corner points" distance from origin
                i[j] = 0.5 + 1 / (2.0 * sqrt(rvs)) * (-1)^(1 + (i[j]))
            end
        end
    end

    #setting axial points
    for i in 1:rvs
        for j in 1:2
            axial_points[2 * i + j - 2, i] = j % 2
        end
    end

    return vcat(two_lvl_points, axial_points)
end
