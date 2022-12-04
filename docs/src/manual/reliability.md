# Reliability Analysis

## Subset Simulation

Subset simulation ([auEstimationSmallFailure2001](@cite)) is an advanced simulation technique for the estimation of small failure probabilities. Here we solve a simple problem where the response $y$ depends on two independent random variables $x_1$ and $x_2$ following a standard normal distribution. The simple linear model is defined by

$y(x_1,x_2) = x_1 + x_2$

with the failure domain

$F = \{(x_1, x_2) : x_1 + x_2 > 9\}.$

The analytical probability of failure can be calculated as

$pf = 1 - \Phi(\frac{9}{\sqrt(2)}) \approx 1 \times 10^{-10}.$

This example is taken from [zuevSubsetSimulationMethod2015](@cite).

In order to solve this, we start by creating the two random variables and group them in a vector `inputs`.

```@example subset
using UncertaintyQuantification, DataFrames # hide
using Random; Random.seed!(8128) # hide
x1 = RandomVariable(Normal(), :x1)
x2 = RandomVariable(Normal(), :x2)
inputs = [x1, x2]
```

Next we define the model as

```@example subset
y = Model(df -> df.x1 + df.x2, :y)
nothing # hide
```

where the first input is our function (which must accept a `DataFrame`) and the second the `Symbol` for the output variable.

To estimate a failure probability we need a performance which is negative if a failure occurs.

```@example subset
function g(df::DataFrame)
    return 9 .- df.y
end
nothing # hide
```

Finally, we create the [`SubSetSimulation`](@ref) object and compute the probability of failure using a standard Gaussian proposal PDF. The value for the target probability of failure at each intermediate level is set to $0.1$ which is generally accepted as the optimal value.

```@example subset
subset = SubSetSimulation(1000, 0.1, 10, Normal())
pf, cov, samples = probability_of_failure(y, g, inputs, subset)

println("Probability of failure: $pf")
```

Alternatively, instead of using the standard Subset simulation algorithm (which internally uses Markov Chain Monte Carlo), we can use [`SubSetInfinity`](@ref) to compute the probability of failure, see [auRareEventSimulation2016](@cite). Here we use a standard deviation of $0.5$ to create the proposal samples for the next level.

```@example subset
subset = SubSetInfinity(1000, 0.1, 10, 0.5)
pf, cov, samples = probability_of_failure(y, g, inputs, subset)

println("Probability of failure: $pf")
```
