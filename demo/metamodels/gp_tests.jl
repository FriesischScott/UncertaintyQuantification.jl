using UncertaintyQuantification
using GaussianProcesses # do we reexport for mean and kernel functions etc.?
using DataFrames
using Random
using Statistics

Random.seed!(20140430)
# Training data
n=10;                          #number of training points
x = 2π * rand(n);              #predictors
y = sin.(x) + 0.05*randn(n) .+ 1e3;   #regressors

#Select mean and covariance function
mZero = MeanZero()                   #Zero mean function
mConst = MeanConst(1.0)
kern = SE(0.0,0.0)                   #Sqaured exponential kernel (note that hyperparameters are on the log scale)

inputs = :x
output = :y
df = DataFrame(inputs => x)
df[!, output] = y

logObsNoise = -1.0                        # log standard deviation of observation noise (this is optional)

gp, = gaussianprocess(df, [inputs], output, kern, mZero, logObsNoise)

# gp = GP(x,y,mConst,kern,logObsNoise)       #Fit the GP
# gp_scaled = GP(x,y_scaled,mZero,kern,logObsNoise-log(std(y))) 

# gp = GP(x,y,mConst,kern)       #Fit the GP
# gp_scaled = GP(x,y_scaled,mZero,kern, -2-log(std(y))) 

# μ, σ² = predict_y(gp,range(0,stop=2π,length=100))
# a, b = predict_y(gp_scaled,range(0,stop=2π,length=100))
# a_ = a .* std(y) .+ mean(y)
# b_ = b .* std(y)

using Optim

optimize_hyperparams!(gp; method=ConjugateGradient(), noise=false)   # Optimise the hyperparameters
# optimize!(gp_scaled; method=ConjugateGradient(), noise=false)

# plot(gp; legend=false, fmt=:png)   #Plot the GP after the hyperparameters have been optimised 

# optimize!(gp; kern = false)   # Don't optimize kernel hyperparameters
# optimize!(gp; kernbounds = [[-1, -1], [1, 1]]) # Optimize the kernel parameters in a box with lower bounds [-1, -1] and upper bounds [1, 1]

# using Plots  #Load Plots.jl package

# scatter(x, y)
# plot!(range(0,stop=2π,length=100), μ, ribbon=σ²)
# plot!(range(0,stop=2π,length=100), a_, ribbon=b_)

# plot(gp; xlabel="x", ylabel="y", title="Gaussian process", legend=false, fmt=:png)      # Plot the GP
# plot(gp_scaled; xlabel="x", ylabel="y", title="Gaussian process", legend=false, fmt=:png)
