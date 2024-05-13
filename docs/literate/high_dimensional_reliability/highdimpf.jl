#===
# High dimensional Subset simulation
## Subset simulation

The implemented subset simulation algorithms `SubSetSimulation`(using component-wise MCMC), `SubSetInfinity` (conditional sampling MCMC) `SubSetInfinityAdaptive` (adaptive conditional sampling MCMC), work efficiently in high dimensions. This tutorial shows how these algorithms scale in dimension `Num_dim` and target probability of failure `pf_target`.

===#

#===
## Example function

In this example, the test model will be sum of independent standard normal distributions

$f_N(X) = \sum^N_i X_i$

$\newline$ where $X_i \sim \Phi(0, 1)$ are standard normal random variables. We will define a linear limitstate: $\newline$ 

$g_N(X) = C_N - f_N(X)$

where $C_N$ will be defined such that the failure probability $\mathbb{P}(g(X) \leq 0)$ matches a pre-defined value `pf_target`. We can find $C_N$ analytically, depending on the `Num_dim` and `pf_target`:

$C_N = F_{\Phi_{\sqrt{N}}}^{-1}(1 - p_{\text{target}})$

where $F_{\Phi_{\sqrt{N}}}^{-1}$ is the quantile function of a Gaussian distribution, with zero mean and variance `sqrt(Num_dim)`.


Since the dimension and failure probability are two parameters of this numerical experiement, we need a compact way to generate an aribitrary number of distributions. This can be done using Julia's metaprogramming.
===#

using UncertaintyQuantification

Num_dim = 2000

for j = 1:Num_dim
    eval(:($(Symbol(:X,j)) = RandomVariable(Normal(0, 1), (Symbol(:X,$j))) ))
end

inputs = [eval(:($(Symbol(:X,i)))) for i = 1:Num_dim]

#===
A model of arbitrary dimensions can be produced using DataFrames, summing the rows corresponding to samples of random variables.
===#

f = Model(
    df -> sum.(eachrow(df[:, names(inputs)])),
    :f,
);

#===
Define the `pf_traget`, C, and the limitstate function
===#

pf_target = 1e-9

fail_limit = quantile(Normal(0, sqrt(Num_dim)), 1 - pf_target)

function g(df)
    return fail_limit .- reduce(vcat, df.f)
end;

#===
Define three simulation methods, and simulate:
===#

simulation_method_1 = SubSetSimulation(2000, 0.1, 20, Uniform(-0.5, 0.5))
simulation_method_2 = SubSetInfinity(2000, 0.1, 20, 0.5)
simulation_method_3 = SubSetInfinityAdaptive(2000, 0.1, 20, 200)

@time pf_1, std_1, samples_1 = probability_of_failure(f, g, inputs, simulation_method_1)
@time pf_2, std_2, samples_2 = probability_of_failure(f, g, inputs, simulation_method_2)
@time pf_3, std_3, samples_3 = probability_of_failure(f, g, inputs, simulation_method_3);

println("True pf: $pf_target | SS: $pf_1 ± $(1.96 * std_1) | SS_inf: $pf_2 ± $(1.96 * std_2) | SS_inf_a: $pf_3 ± $(1.96 * std_3)")

#===
We note that although basic Monte Carlo simulation works independently of dimension, for a target failure probability of $10^{-9}$, even with a billion $10^9$ samples can give $p_{\text{f}}=0$.
===#



function run_sim(Num_dim, pf_target, N_runs, N_batches = 50) #hide

    for j = 1:Num_dim                                                           #hide
        eval(:($(Symbol(:X,j)) = RandomVariable(Normal(0, 1), (Symbol(:X,$j))) ))#hide
    end#hide

    inputs = [eval(:($(Symbol(:X,i)))) for i = 1:Num_dim]                       #hide

    y = Model(                                                                  #hide
        df -> sum.(eachrow(df[:, names(inputs)])),#hide
        :y,#hide
    )#hide

    fail_limit = quantile(Normal(0, sqrt(Num_dim)), 1 - pf_target)#hide

    function limitstate(df)#hide
        return fail_limit .- reduce(vcat, df.y)#hide
    end#hide

    simulation_method_1 = SubSetSimulation(N_runs, 0.1, 20, Uniform(-0.5, 0.5))#hide
    simulation_method_2 = SubSetInfinity(N_runs, 0.1, 20, 0.5)#hide
    simulation_method_3 = SubSetInfinityAdaptive(N_runs, 0.1, 20, Integer(floor(N_runs * 0.1)))#hide

    pf_1 = zeros(N_batches)#hide
    pf_2 = zeros(N_batches)#hide
    pf_3 = zeros(N_batches)#hide

    std_1 = zeros(N_batches)#hide
    std_2 = zeros(N_batches)#hide
    std_3 = zeros(N_batches)#hide

    for i = 1:N_batches#hide
        pf_1[i], std_1[i], _ = probability_of_failure(y, limitstate, inputs, simulation_method_1)#hide
        pf_2[i], std_2[i], _ = probability_of_failure(y, limitstate, inputs, simulation_method_2)#hide
        pf_3[i], std_3[i], _ = probability_of_failure(y, limitstate, inputs, simulation_method_3)#hide
    end#hide

    return pf_1, pf_2, pf_3, std_1, std_2, std_3#hide
end;#hide

using Plots#hide

pfs = 10.0 .^collect(-1:-1:-13)#hide
N_dimension = 200#hide
N_samples = 2000#hide
N_batchs = 1#hide

pf_1 = zeros(length(pfs), N_batchs)#hide
pf_2 = zeros(length(pfs), N_batchs)#hide
pf_3 = zeros(length(pfs), N_batchs)#hide
#hide
for (i, pf) in enumerate(pfs)#hide
    (pf_1[i,:], pf_2[i,:], pf_3[i,:], _, _, _) = run_sim(N_dimension, pf, N_samples, N_batchs)#hide
end#hide
#hide
pf_1_mean = mean(pf_1, dims = 2)#hide
pf_2_mean = mean(pf_2, dims = 2)#hide
pf_3_mean = mean(pf_3, dims = 2)#hide
#hide
scaleit(x) = -log10.(x)#hide


plot(scaleit(pfs), scaleit(pf_1_mean),  seriestype=:scatter, label="SubSet")#hide
plot!(scaleit(pfs), scaleit(pf_2_mean), seriestype=:scatter, label="SubSetInfinity")#hide
plot!(scaleit(pfs), scaleit(pf_3_mean), seriestype=:scatter, label="SubSetAdaptive", legend=:topleft, xlabel ="-log10(pf) target",  ylabel ="-log10(pf) simulated")#hide
plot!(scaleit(pfs), scaleit(pfs),  minorgrid=true, label = false, lw = 2)#hide
title!("Num_dims = $N_dimension, N_samples = $N_samples")#hide

pf = 10^-4#hide
N_dimension = 200#hide
N_samples = Integer.(floor.(2 .^(range(1,stop=10,length=10))*10))#hide
N_batchs = 20#hide

pf_1_samples = zeros(length(N_samples), N_batchs)#hide
pf_2_samples = zeros(length(N_samples), N_batchs)#hide
pf_3_samples = zeros(length(N_samples), N_batchs)#hide

for (i, N_s) in enumerate(N_samples)#hide
    (pf_1_samples[i,:], pf_2_samples[i,:], pf_3_samples[i,:], _, _, _) = run_sim(N_dimension, pf, N_s, N_batchs)#hide
end#hide

pf_1_samples_mean = mean(pf_1_samples, dims = 2)#hide
pf_2_samples_mean = mean(pf_2_samples, dims = 2)#hide
pf_3_samples_mean = mean(pf_3_samples, dims = 2)#hide

pf_1_samples_lo = [quantile(pf_1_samples[i, :], 0.025) for i = 1:size(pf_1_samples, 1)]#hide
pf_1_samples_hi = [quantile(pf_1_samples[i, :], 0.975) for i = 1:size(pf_1_samples, 1)]#hide

pf_2_samples_lo = [quantile(pf_2_samples[i, :], 0.025) for i = 1:size(pf_2_samples, 1)]#hide
pf_2_samples_hi = [quantile(pf_2_samples[i, :], 0.975) for i = 1:size(pf_2_samples, 1)]#hide

pf_3_samples_lo = [quantile(pf_3_samples[i, :], 0.025) for i = 1:size(pf_3_samples, 1)]#hide
pf_3_samples_hi = [quantile(pf_3_samples[i, :], 0.975) for i = 1:size(pf_3_samples, 1)]#hide


plot(N_samples, pf_1_samples_hi, fillrange = pf_1_samples_lo, label="SubSet", alpha = 0.2, color = theme_palette(:auto)[1])#hide
plot!(N_samples, pf_2_samples_hi, fillrange = pf_2_samples_lo, label="SubSetInfinity", alpha = 0.2, color = theme_palette(:auto)[2])#hide
plot!(N_samples, pf_3_samples_hi, fillrange = pf_3_samples_lo, label="SubSetAdaptive", alpha = 0.2, color = theme_palette(:auto)[3])#hide

plot!(N_samples, pf_1_samples_mean, color = theme_palette(:auto)[1], label = false, lw = 2, seriestype=:scatter)#hide
plot!(N_samples, pf_2_samples_mean, color = theme_palette(:auto)[2], label = false, lw = 2, seriestype=:scatter)#hide
plot!(N_samples, pf_3_samples_mean, color = theme_palette(:auto)[3], label = false, lw = 2, seriestype=:scatter)#hide

plot!(N_samples, fill(pf, length(N_samples)),  minorgrid=true, xscale=:log10, yscale=:log10, label="target pf", legend=:bottomright, xlabel ="samples per level", ylabel="pf", lw = 2, color = theme_palette(:auto)[4])#hide
title!("Num_dims = $N_dimension, targe pf = $pf")#hide


pf = 10^-4#hide
N_dimension = Integer.(floor.(2 .^(range(1,stop=10,length=10))*10))#hide
N_samples = 2000#hide
N_batchs = 1#hide

pf_1_dims = zeros(length(N_dimension), N_batchs)#hide
pf_2_dims = zeros(length(N_dimension), N_batchs)#hide
pf_3_dims = zeros(length(N_dimension), N_batchs)#hide

std_1_dims = zeros(length(N_dimension), N_batchs)#hide
std_2_dims = zeros(length(N_dimension), N_batchs)#hide
std_3_dims = zeros(length(N_dimension), N_batchs)#hide

for (i, Ndim) in enumerate(N_dimension)#hide
    (pf_1_dims[i,:], pf_2_dims[i,:], pf_3_dims[i,:], std_1_dims[i,:], std_2_dims[i,:], std_3_dims[i,:]) = run_sim(Ndim, pf, N_samples, N_batchs)#hide
end#hide


pf_1_dims_mean = mean(pf_1_dims, dims = 2)#hide
pf_2_dims_mean = mean(pf_2_dims, dims = 2)#hide
pf_3_dims_mean = mean(pf_3_dims, dims = 2)#hide

std_1_dims_mean = mean(std_1_dims, dims = 2)#hide
std_2_dims_mean = mean(std_2_dims, dims = 2)#hide
std_3_dims_mean = mean(std_3_dims, dims = 2)#hide


plot(N_dimension, pf_1_dims_mean - 1.96 * std_1_dims_mean, fillrange = pf_1_dims_mean + 1.96 * std_1_dims_mean, label="SubSet", alpha = 0.2, color = theme_palette(:auto)[1])#hide
plot!(N_dimension, pf_2_dims_mean - 1.96 * std_2_dims_mean, fillrange = pf_2_dims_mean + 1.96 * std_2_dims_mean, label="SubSetInfinity", alpha = 0.2, color = theme_palette(:auto)[2])#hide
plot!(N_dimension, pf_3_dims_mean - 1.96 * std_3_dims_mean, fillrange = pf_3_dims_mean + 1.96 * std_3_dims_mean, label="SubSetAdaptive", alpha = 0.2, color = theme_palette(:auto)[3])#hide

plot!(N_dimension, pf_1_dims_mean,  seriestype=:scatter, color = theme_palette(:auto)[1], label = false)#hide
plot!(N_dimension, pf_2_dims_mean, seriestype=:scatter, color = theme_palette(:auto)[2], label = false)#hide
plot!(N_dimension, pf_3_dims_mean, seriestype=:scatter, color = theme_palette(:auto)[3], label = false)#hide
plot!(N_dimension, fill(pf, 10),  minorgrid=true, label = false, xscale=:log10, legend=:topleft, xlabel ="N dimensions",  ylabel ="pf")#hide
title!("N_samples = $N_samples, pf target = $pf")#hide