using UncertaintyQuantification
using DataFrames
using KernelDensity
using Plots
x = RandomVariable(Normal(), :x)

samples = sample(x, 100) # sample(x, MonteCarlo(100))

bin = Array({Int64})

"User defined bin size"
"Warning when bin beyond 50"
###Biggest bin size
##Adding bin parameter
function binning(bin_size::Int64, to_bin::DataFrame)
    if (nrow(to_bin)% bin_size == 0)
        times = nrow(to_bin) ÷ bin_size
        bin_num = 0
        bin =  Vector{Int64}()
        for i in 1:times
        append!(bin,repeat([bin_num], outer = bin_size))
        bin_num += 1
        end
    else
        counter = 0
        bin_num = 1
        bin =  Vector{Int64}()
        while counter != nrow(to_bin)
            for i in 1:bin_size
            push!(bin,bin_num)
            end
            bin_num+=1
            counter += bin_size
            println(counter)   
        end
    end
return bin
end


function sequential_binning(bin_size::Int64, to_bin::DataFrame)
    times = nrow(to_bin) ÷ bin_size + 1
    bin =  Vector{Int64}()
    for i in 1:times
        append!(bin,repeat(collect(1:1:nrow(to_bin) ÷ bin_size), outer = bin_size))
    end
return bin[1:nrow(to_bin)]
end


function sk_extraction(df::DataFrame)
    vec = []
    for i in 1:nrow(unique(samples, :bin))
        samples[(samples.bin .= 1), :]
        push!(vec,kde(df[(df.bin .= i), :].x))
    end
return vec
end



db_cor_estim(s, n)

end