using LinearAlgebra, DataFrames, FiniteDifferences, Dierckx, Reexport,
Distributions, Random, Plots

mcit = 10000

new_x(x) = Normal{Float64}(x[2]/2, sqrt(3/4))
new_y(x) = Normal{Float64}(x[1]/2, sqrt(3/4))

condist = [new_x, new_y]
ample = zeros(mcit,size(condist,1))
ample[1,:] = [1 -1]
dsample = DataFrame(x = ample[:,1], y = ample[:,2])

function gibbssample(
    condist::Array,
    startval::Array,
    mcit::Int64,
    burn::Int64
    )

    sample = zeros(mcit, size(startval, 2))
    sample[1,:] = startval

    for i in 1:mcit-1
        sample[i+1,:] = sample[i,:]
        for j in 1:size(sample, 2)
            sample[i+1,j] = rand(condist[j](sample[i+1,:]))
        end
    end
    return sample[burn+1:mcit,:]
end

dsample = gibbssample(condist, [1 -1], mcit, 500)

display(scatter(dsample[:,1], dsample[:, 2]))
#display(histogram(dsample[:, 1]))


#calc mean
function samplemean(
    dsample::Array
    )
    sumx = sum(dsample, dims=1)
    mcit = size(dsample,1)
    return sumx./mcit
end

mean1 = samplemean(dsample)


#calc cedible interval



function calcequaltails(
    dsample,
    quantile::Float64,
    )
    snum = size(dsample, 1)
    dsample = sort(dsample, dims = 1)

    low = floor(Int32,snum*(1-quantile)/2)
    high = floor(Int32,snum - snum*(1-quantile)/2)

    rv = zeros(2, size(dsample,2))

    for i in 1:size(dsample,2)
        rv[:,i] = [dsample[low,i]; dsample[high,i]]
    end
    return rv
end

eqt = calcequaltails(dsample, 0.95)
