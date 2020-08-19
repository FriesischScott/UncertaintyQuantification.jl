function sobolindices(
    models::Array{<:UQModel},
    input::Array{<:UQInput},
    sim::MonteCarlo
)
    # the final model that is to be evaluated has to be named :fkt

    asamples = sample(input, sim.n)
    bsamples = sample(input, sim.n)

    numberofinputs = 0
    for i in 1:size(asamples, 2)
        if isa(input[i], RandomVariable)
            numberofinputs += 1
        end
    end


    for m in models
        evaluate!(m, asamples) # calculation of the model using the random values in matrix A
        evaluate!(m, bsamples) # calculation of the model using the random values in matrix B
    end

    V_i = Vector{Float64}(undef, numberofinputs)
    randvarcounter = 0

    for i in 1:length(input)
        if isa(input[i], RandomVariable)    # only randomvariable are swapped
            abisamples = DataFrame()
            randvarcounter += 1

            for j in 1:length(input)
                if j == i
                    abisamples[!, input[j].name] = bsamples[!,j] # each random variable column is swapped once with a bsamples column
                else
                    abisamples[!, input[j].name] = asamples[!,j]
                end
            end

            for m in models
                evaluate!(m, abisamples)  # recalculating the model with the new column
            end

            V_i[randvarcounter] = 1/(2*sim.n) * (sum((asamples.fkt .- abisamples.fkt).^2)) # estimator (f) from Saltelli et al
        end
    end

    f_0 = 1/sim.n * sum((asamples.fkt .* bsamples.fkt)) # given estimator for (f_0)^2

    V_Y = 1/sim.n * sum((asamples.fkt).^2) - f_0  # given estimator for V(Y)
    out = V_i ./ V_Y # normalize the values using V(Y)
    return out
end
