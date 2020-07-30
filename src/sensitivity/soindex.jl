function soindex(
    models::Array{<:UQModel},
    input::Array{<:UQInput},
    sim::MonteCarlo
)
    # each variable input i needs two entries named Ai, Bi
    # other input names are free to choose as long as they are correctly called upon by the model
    # the final model that is to be evaluated has to be named :fkt

    samples = sample(input, sim.n)
    msamples = copy(samples)

    numberofinputs = convert(Int32, size(samples, 2)/2)

    for m in models
        evaluate!(m, samples) # calculation of the model using the random values in matrix A
    end

    V_i = Vector{Float64}(undef, numberofinputs)

    for i in 1:numberofinputs
        rename!(msamples, "A$i" => "X")
        rename!(msamples, "B$i"=>"A$i") # changing the df name so that the model now uses the i-th column of matrix B

        for m in models
            evaluate!(m, msamples)  # recalculating the model with the new column
        end

        V_i[i] = 1/(2*sim.n) * (sum((samples.fkt .- msamples.fkt).^2)) # estimator (f) from Saltelli et al

        rename!(msamples, "A$i" => "OUT$i") # this column is no longer used
        rename!(msamples, "X"=>"A$i")   # this column will still be used for the next iteration
    end

    msamples = copy(samples)

    for i in 1:numberofinputs
        rename!(msamples, "A$i" => "OUT$i") # swapping all columns of A with B for (f_0)^2 estimator
        rename!(msamples, "B$i"=>"A$i")
    end

    for m in models
        evaluate!(m, msamples) # recalculating the model with matrix B
    end

    f_0 = 1/sim.n * sum((samples.fkt .* msamples.fkt)) # given estimator for (f_0)^2

    V_Y = 1/sim.n * sum((samples.fkt).^2) - f_0  # given estimator for V(Y)
    out = V_i ./ V_Y # normalize the values using V(Y)
    return out
end
