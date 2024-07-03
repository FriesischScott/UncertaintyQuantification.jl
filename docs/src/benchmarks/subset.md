# High dimensional Subset simulation

## Subset simulation

The implemented subset simulation algorithms [`SubSetSimulation`](@ref)(using component-wise MCMC), [`SubSetInfinity`](@ref)(conditional sampling MCMC) and [`SubSetInfinityAdaptive`](@ref)(adaptive conditional sampling MCMC), work efficiently in high dimensions. This benchmark shows how these algorithms scale with increasing number of dimension `N` and increasingly smaller target probability of failure `pf_target`.

## Example function

In this example, the test model will be sum of independent standard normal distributions

```math
f_N(X) = \sum^N_i X_i,
```

where $X_i \sim \Phi(0, 1)$ are standard normal random variables. We will define a linear limitstate

```math
g_N(X) = C_N - f_N(X),
```

where $C_N$ will be defined such that the failure probability $\mathbb{P}(g(X) \leq 0)$ matches a pre-defined value `pf_target`.

We can find $C_N$ analytically, depending on the chosen number of dimensions and target probability of failure

```math
C_N = F_{\Phi_{\sqrt{N}}}^{-1}(1 - pf_{\text{target}}),
```

where $F_{\Phi_{\sqrt{N}}}^{-1}$ is the quantile function of a Gaussian distribution, with zero mean and standard deviation `sqrt(N)`.

Since the dimension and failure probability are two parameters of this numerical experiment, we can dynamically create the required number of random variables using broadcasting

```julia
using UncertaintyQuantification

N = 2000

 inputs = RandomVariable.(Normal(), [Symbol("x$i") for i in 1:N])

```

The model can be defined generalized for arbitrary dimensions by summing the columns of the `DataFrame`. Using `names(inputs)` to select the columns we can safely exclude any extra variables that might be present.

```julia

f = Model(
    df -> sum(eachcol(df[:, names(inputs)])),
    :f
)
```

Next, the `pf_target` and corresponding limit state are defined.

```julia

pf_target = 1e-9

fail_limit = quantile(Normal(0, sqrt(N)), 1 - pf_target)

function g(df)
    return fail_limit .- df.f
end
```

For this benchmark, the probability of failure will be estimated using all available variants of Subset simulation

```julia
simulation_method_1 = SubSetSimulation(2000, 0.1, 20, Uniform(-0.5, 0.5))
simulation_method_2 = SubSetInfinity(2000, 0.1, 20, 0.5)
simulation_method_3 = SubSetInfinityAdaptive(2000, 0.1, 20, 200)

pf_1, std_1, samples_1 = probability_of_failure(f, g, inputs, simulation_method_1)
pf_2, std_2, samples_2 = probability_of_failure(f, g, inputs, simulation_method_2)
pf_3, std_3, samples_3 = probability_of_failure(f, g, inputs, simulation_method_3)

println("True pf: $pf_target | SS: $pf_1 ± $(1.96 * std_1) | SS_inf: $pf_2 ± $(1.96 * std_2) | SS_inf_a: $pf_3 ± $(1.96 * std_3)")
```

!!! note "Monte Carlo simulation"
    Although standard Monte Carlo simulation works independently of dimension, for a target failure probability of $10^{-9}$, even with a billion $10^9$ samples it can return $p_f=0$.
