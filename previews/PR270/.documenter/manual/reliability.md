
# Reliability Analysis {#Reliability-Analysis}

In the context of structural engineering, engineering design, and risk assessment, the term reliability is used to describe the ability of system to perform its intended function under varying conditions over time.

There, the state of a system is identified by its _performance function_ (also sometimes refered to as limit state function) $g(\boldsymbol{x})$ such that:

$$g(\boldsymbol{x}) =
\begin{cases}
    > 0 & \text{safe\ domain}\\
    \leq 0 & \text{failure \ domain}\\
\end{cases}.$$

Then the _probability of failure_ is defined as the likelihood of the system being in the failed state, given as

$$p_f = \int_{g(\boldsymbol{x}) \leq 0} f_{\boldsymbol{X}}(\boldsymbol{x}) \mathrm{d} \boldsymbol{x}.$$

Here, $f_{\boldsymbol{X}}(\boldsymbol{x})$ denotes the joint probability density function (PDF) of the input $\boldsymbol{X}$.

## Definition of the Input, Model and Performance {#Definition-of-the-Input,-Model-and-Performance}

The first step of the implementation of a reliability analysis in `UncertaintyQuantification.jl` is the definition of the probabilistic input and the model which is shown exemplarily.

Here we use a modification of the first example presented in [[12](/references#papaioannouCombination2021)] which uses a quadratic performance function and a probabilistic input containing of two standard normal random variables $\boldsymbol{X} = [X_1, X_2]$ with $X_i \sim \mathcal{N}(0,1)$. The model is given as

$$y(\boldsymbol{X}) = 0.1(X_1 - X_2)^2 - \frac{1}{\sqrt{2}} (X_1 + X_2).$$

Then, the performance function is defined as

$$g(\boldsymbol{X}) = y(\boldsymbol{X}) + 4.$$

The probabilistic input is implemented as

```julia
x = RandomVariable.(Normal(), [:x1, :x2])
```


```ansi
2-element Vector{RandomVariable{Normal{Float64}}}:
 RandomVariable{Normal{Float64}}(Normal{Float64}(μ=0.0, σ=1.0), :x1)
 RandomVariable{Normal{Float64}}(Normal{Float64}(μ=0.0, σ=1.0), :x2)
```


Next we define the model for the response $y(\boldsymbol{X})$ as

```julia
y = Model(df -> 0.1*(df.x1 - df.x2).^2 - 1/sqrt(2) * (df.x1 + df.x2), :y)
```


where the first input is the function $y$ (which must accept a `DataFrame`) and the second argument is the `Symbol` for the output variable. With the help of the model, we can define the performance function $g$ which again takes a `DataFrame` as an input:

```julia
g(df) = df.y .+ 4
```


## Approximation Methods {#Approximation-Methods}

### First Order Reliability Method {#First-Order-Reliability-Method}

The First Order Reliability Method (FORM) [[13](/references#rackwitzStructuralReliability1978)] estimates the failure probability by finding a linear approximation of the performance function at the so-called _design point_ $\boldsymbol{U}^*$. The design point represents the point on the surface of the performance function $g(\boldsymbol{X}) = 0$ that is closest to the origin in the standard normal space.

That distance from the design point to the origin is referred to as the _reliability index_ given as $\beta^* = ||\boldsymbol{U}^*||$. Due to the transformation to the standard normal space, the probability of failure is simply given as

$$\hat{p}_{f, \mathrm{FORM}} = \Phi(-\beta^*)$$

where $\Phi$ denotes the standard normal CDF.

In addition to the $\beta^*$, the location of the design point is specified by the _important direction_ defined as:

$$\boldsymbol{\alpha}^* = \frac{\boldsymbol{U}^*}{||\boldsymbol{U}^*||}.$$

In `UncertaintyQuantification.jl` a FORM analysis can be performed calling `probability_of_failure(model, performance, input, simulation)` where `FORM()` is passed as the simulation method:

```julia
pf_form, β, dp = probability_of_failure(y, g, x, FORM())

println("Probability of failure: $pf_form")
```


```ansi
Probability of failure: 3.1671241833194436e-5
```


## Simulation Methods {#Simulation-Methods}

### Monte Carlo Simulation {#Monte-Carlo-Simulation}

Monte Carlo Simulation (MCS) offers an approximation of the failure probability using stochastic simulation.

It utilizes an indicator function of the failure domain

$$\mathbb{I}[g(\boldsymbol{x})] =
\begin{cases}
    0 & \text{when} \ g(\boldsymbol{x}) > 0\\
    1 & \text{when} \ g(\boldsymbol{x}) \leq 0\\
\end{cases}.$$

This allows for the failure probability to be interpreted as the expected value of the indicator function

$$p_f = \int_{\boldsymbol{X}} \mathbb{I}[g(\boldsymbol{x})]\ f_{\boldsymbol{X}}(\boldsymbol{x}) \mathrm{d} \boldsymbol{x} = \mathbb{E}\big[\mathbb{I}[g(\boldsymbol{x})]\big].$$

The Monte Carlo estimate of the failure probability is given as

$$p_f \approx \hat{p}_f = \frac{1}{N} \sum_{i=1}^N \mathbb{I}[g(\boldsymbol{x}_i)]$$

where $\{\boldsymbol{x}_i\}_{i=1}^N$ represents a set of $N$ samples drawn from the input PDF $f_{\boldsymbol{X}}(\boldsymbol{x})$. The variance of the estimator is given as

$$\operatorname{Var}[\hat{p}_f] = \frac{\hat{p}_f (1-\hat{p}_f)}{N}.$$

In `UncertaintyQuantification.jl` we can perform a Monte Carlo Simulation by defining the analysis as `MonteCarlo(n)`  where `n` is the number of samples:

```julia
mc = MonteCarlo(10^7)
```


Then the reliability analysis is performed by calling `probability_of_failure(model, performance, input, simulation)`.

```julia
pf_mc, std_mc, samples = probability_of_failure(y, g, x, mc)

println("Probability of failure: $pf_mc")
println("Coefficient of variation: $(std_mc/pf_mc)")
```


```ansi
Probability of failure: 1.82e-5
Coefficient of variation: 0.07412425712616279
```


### Importance Sampling {#Importance-Sampling}

Based on the standard MCS method, a class of advanced method exist that have to goal to accelerate the estimation of the failure probability by requiring fewer model calls.

Importance Sampling [[14](/references#melchersImportanceSampling1989)] introduces a second density that is _biased_ in a way that it generates more samples in the failure domain. Typically such a density is constructed around the design point obtained in a preceding FORM analysis.

In order to perform a reliability analysis using Importance Sampling, we again have to specify the number of samples and then can `probability_of_failure()`.

```julia
is = ImportanceSampling(1000)
pf_is, std_is, samples = probability_of_failure(y, g, x, is)

println("Probability of failure: $pf_is")
println("Coefficient of variation: $(std_is/pf_is)")
```


```ansi
Probability of failure: 2.0801590194268694e-5
Coefficient of variation: 0.04739109089805725
```


### Radial Based Importance Sampling {#Radial-Based-Importance-Sampling}

Radial based importance sampling (RBIS) [[15](/references#harbitzEfficientSamplingMethod1986)] increases the efficiency of the Monte Carlo simulation by sampling in standard normal space and excluding a β-sphere where no failures occur from the sampling domain. Here, `β` is the reliability index obtained from a preliminary analysis like FORM. The probability of failure is then estimated as

$$p_f \approx \hat{p}_f = (1- \chi^2_k(\beta^2) \frac{1}{N} \sum_{i=1}^N \mathbb{I}[g(\boldsymbol{x}_i)],$$

where $\chi^2_k$ is the CDF of the Chi-squared distribution with $k$ degrees of freedom and $k$ is the number of random variables.

If no `β` or `β=0.0` is passed to the [`RadialBasedImportanceSampling`](/api/simulations#UncertaintyQuantification.RadialBasedImportanceSampling) constructor, a FORM analysis will automatically be performed.

```julia
rbis = RadialBasedImportanceSampling(1000)
pf_rbis, std_rbis, samples = probability_of_failure(y, g, x, rbis)

println("Probability of failure: $pf_rbis")
println("Coefficient of variation: $(std_rbis/pf_rbis)")
```


```ansi
Probability of failure: 1.9456832418391024e-5
Coefficient of variation: 0.12744167022738218
```


A scatter plot clearly shows the exclusion of the β-sphere. 
![](rbis-samples.svg)


### Line Sampling {#Line-Sampling}

Another advanced Monte Carlo method for reliability analysis is Line Sampling [[16](/references#koutsourelakisReliability2004)]. Its main idea is to use parallel lines for sampling rather than points.

Therefore first the problem is transformed into the standard normal space to make use of the invariance of rotation. The important direction $\boldsymbol{\alpha}$ is determined, e.g., using FORM or the gradient at the origin. Then, samples are generated and projected onto the hyperplane orthogonal to $\boldsymbol{\alpha}$. From each point on the hyperplane, a line is drawn parallel to $\boldsymbol{\alpha}$ and its intersection with the performance function is determined using root finding based on a spline interpolation scheme, giving the set of distances $\{\beta^{(i)}\}_{i=1}^N$ from the hyperplane to the intersection with the performance function. Due to working in the standard normal space, the _failure probability along each line_ is given as

$$p_{f, \mathrm{line}}^{(i)} = \Phi(-\beta^{(i)})$$

Finally, the probability of failure is obtained as the mean of the failure probabilities along the lines

$$\hat{p}_{f,\mathrm{LS}} = \frac{1}{N} \sum_{i=1}^N p_{f, \mathrm{line}}^{(i)}.$$

The variance of $\hat{p}_{f,\mathrm{LS}}$ is given by the variance of the line failure probabilities:

$$\operatorname{Var}[\hat{p}_{f,\mathrm{LS}}] = \frac{1}{N(N-1)} \sum_{i=1}^N \Big(p_{f, \mathrm{line}}^{(i)} - \hat{p}_{f,\mathrm{LS}}\Big)^2.$$

Similar to standard MCS, we have to pass $N$ to the Line Sampling method. However, here we pass the number of lines. Optionally, we can pass a vector of the points along each line that are used to evaluate the performance function and a predetermined direction $\boldsymbol{\alpha}$:

```julia
ls = LineSampling(100, collect(0.5:0.5:8.0))
pf_ls, std_ls, samples = probability_of_failure([y], g, x, ls)

println("Probability of failure: $pf_ls")
println("Coefficient of variation: $(std_ls/pf_ls)")
```


```ansi
Probability of failure: 1.7598929131182565e-5
Coefficient of variation: 0.05972209531902246
```


### Advanced Line Sampling {#Advanced-Line-Sampling}

Advanced Line Sampling [[17](/references#deangelisAdvances2015)] is a further enhancement of the standard line sampling methods due to two main features:
1. The important direction $\boldsymbol{\alpha}$ is adapted once a more probable point is found
  
2. The lines are processed sorted by proximity of the points on the hyperplane.
  

Especially the second point enables the use of an iterative root finder using Newton&#39;s method.

The definition of the `AdvancedLineSampling` simulation method is similar to that of regular Line Sampling. The number of lines has to be given to the constructor and we can optionally give the number of points along the line which is only used to find the starting point of the iterative root search.

```julia
als = AdvancedLineSampling(100, collect(0.5:0.5:10))
pf_als, std_als, samples = probability_of_failure([y], g, x, als)

println("Probability of failure: $pf_als")
println("Coefficient of variation: $(std_als/pf_als)")
```


```ansi
Probability of failure: 1.923919866421047e-5
Coefficient of variation: 0.05401467892474384
```


For `AdvancedLineSampling`, we can also define the (initial) direction and options of the iterative root finding, i.e., the `tolerance`, `stepsize` of the gradient and `maxiterations`.

::: tip Parallelism

We note that Advanced Line Sampling is a serial algorithm, although much fewer samples (order of magnitude) are required. If a large amount of parallel compute is available, standard Line Sampling may be more attractive, which is &quot;embarrassingly&quot; parallel like Monte Carlo.

:::

## Subset Simulation {#Subset-Simulation}

Subset simulation [[18](/references#auEstimationSmallFailure2001)] is an advanced simulation technique for the estimation of small failure probabilities. This approach involves decomposing the problem into a sequence of conditional probabilities that are estimated using Markov Chain Monte Carlo.

We create the [`SubSetSimulation`](/api/simulations#UncertaintyQuantification.SubSetSimulation) object and compute the probability of failure using a standard Gaussian proposal PDF. The value for the target probability of failure at each intermediate level is set to $0.1$ which is generally accepted as the optimal value.

```julia
subset = SubSetSimulation(1000, 0.1, 10, Normal())
pf_sus, std_sus, samples = probability_of_failure(y, g, x, subset)

println("Probability of failure: $pf_sus")
println("Coefficient of variation: $(std_sus/pf_sus)")
```


```ansi
Probability of failure: 1.4892e-5
Coefficient of variation: 0.34462771589639685
```


Alternatively, instead of using the standard Subset simulation algorithm (which internally uses Markov Chain Monte Carlo), we can use [`SubSetInfinity`](/api/simulations#UncertaintyQuantification.SubSetInfinity) to compute the probability of failure, see [[19](/references#auRareEventSimulation2016)]. Here we use a standard deviation of $0.5$ to create the proposal samples for the next level.

```julia
subset = SubSetInfinity(1000, 0.1, 10, 0.5)
pf_sus, std_sus, samples = probability_of_failure(y, g, x, subset)

println("Probability of failure: $pf_sus")
println("Coefficient of variation: $(std_sus/pf_sus)")
```


```ansi
Probability of failure: 1.4200000000000003e-5
Coefficient of variation: 0.3167927200155807
```


## Imprecise Reliability Analysis {#Imprecise-Reliability-Analysis}

if _epistemic_ uncertainty is considered in the input variables it must also be propagated through the analysis, converting the probability of failure into an interval itself. The goal of the reliability analysis is then to find the lower bound $\underline{p}_f$ and upper bound $\overline{p}_f$ of the probability of failure such that $p_f \in [\underline{p}*f, \overline{p}*f]$. Formally, this can be defined as

$$    \underline{p}*f = \min*{\vartheta \in \mathbf{\theta}} \int \mathbb{I}[g(x)]f*{\mathbf{X}}(x, \vartheta) dx$$

and

$$\overline{p}*f = \max*{\vartheta \in \mathbf{\theta}} \int  \mathbb{I}[g(x)]f*{\mathbf{X}}(x, \vartheta) dx,$$

where $\theta$ represents the epistemic parameters characterizing the inputs.

The following two sections provide two possible solutions to this challenging problem.

### Double-loop Monte Carlo Simulation {#Double-loop-Monte-Carlo-Simulation}

The simplest way to solve the imprecise reliability analysis is by double-loop simulation. The name refers to the need for two loops in comparison to a standard reliability analysis. The outer-loop essentially solves two optimisation problems over the parameter space of the epistemic inputs to find the combinations that minimize and maximize the probability of failure. The inner-loop requires a reliability calculation, with the current combination of parameters under consideration being mapped to precise inputs. In practice, [`IntervalVariable`](/api/inputs#UncertaintyQuantification.IntervalVariable) is mapped to a [`Parameter`](/api/inputs#UncertaintyQuantification.Parameter) while a [`ProbabilityBox`](/api/inputs#UncertaintyQuantification.ProbabilityBox) nested inside a [`RandomVariable`](/api/inputs#UncertaintyQuantification.RandomVariable) yields a regular `RandomVariable`. Therefore, the double-loop simulation treats p-boxes as parametric. Then, a comprehensive reliability analysis using these purely _aleatory_ inputs is carried out. This repeated analysis in the inner-loop makes the double-loop simulation computationally demanding. If a Monte Carlo simulation is applied in the inner-loop this is known as double-loop Monte Carlo simulation.

Special attention must be paid to the type of optimisation algorithm used, as random sampling in the inner-loop leads to non-smooth objective functions. Here we have chosen to apply mesh adaptive direct search (MADS) algorithms which are specifically designed for such cases [[20](/references#abramsonOrthoMADSDeterministicMADS2009)]. However, exploration of alternative methods is part of the ongoing development.

Estimating the probability of failure is effectively separated into two independent problems, one for each bound. This provides the ability to apply different types of analyses. For example, using a larger number of samples for the lower bound.

As an example, consider the function

$$     f(x_1,x_2) = x_1 +x_2,$$

where $x_1 \in [-1,1]$ is an interval and $x_2 \sim N(0, [1,2])$ is distributed as a Gaussian p-box. The associated performance function is

$$     g(x) = 9 + f(x),$$

i.e., failures are $f(x) \leq -9$. The analytical solution to this problem is known to be

$$    \underline{p}_f = \Phi(-9; 1, 1)$$

and

$$     \overline{p}_f = \Phi(-9; -1, 2),$$

where $\Phi(x; \mu, \sigma)$ is the Gaussian CDF with mean $\mu$ and standard deviation $\sigma$. The follow code shows how to set up this problem in _UncertaintyQuantification.jl_.

```julia
     x1 = IntervalVariable(-1,1,:x1)
     x2 = RandomVariable(ProbabilityBox{Normal}(Dict(:μ => 0, :σ => Interval(1,2,))), :x2)
     f = Model(df -> df.x1 .+ df.x2, :f)
```


Then, the relibility analysis is again run through the  [`probability_of_failure`](/api/reliability#UncertaintyQuantification.probability_of_failure-Tuple{Union{UQModel,%20Vector{<:UQModel}},%20Function,%20Union{UQInput,%20Vector{<:UQInput}},%20FORM}) function. As the final argument we we pass the inner simulation wrapped inside a [`DoubleLoop`](/api/reliability#UncertaintyQuantification.DoubleLoop). For this simple example we apply the first order reliability method (FORM), as this can reliably estimate the small failure probability of the lower bound which is (\approx 7.6e-24). The bounds on the `p_\text{f}` can be obtained from the output interval through `pf.lb` and `pf.ub`. The epistemic uncertainty is propagated correctly and matches the analytical solution. On top of the probability of failure, the double-loop analysis also returns the parameters that lead to the lower and upper bound. Here, they are correctly identified as $([1,1]$ and $[-1,2]$.

```julia
     pf, x_lb, x_ub = probability_of_failure(f, df -> 9 .+ df.f, [x1,x2], DoubleLoop(FORM()))
```


```ansi
([7.619853024067617e-24, 3.167124183314543e-5], [1.0, 1.0], [-1.0, 2.0])
```


### Random Slicing {#Random-Slicing}

An alternative method for computing probability bounds for reliability problems is based on random-set theory, as outlined by [[21](/references#alvarez2018estimation)]. We colloquially call this &quot;random slicing&quot;, as will become apparent. As opposed to double-loop Monte Carlo, random slicing is a _distribution-free_ or a _non-parametric_ technique, as it does not make use of distribution parameters or family. For this reason, it is slightly more general, but can provide wider probability intervals in certain simulations.

In random slicing, we make use of the fact that a p-box is a random-set [[5](/references#fersonConstructingProbabilityBoxes2015)] to simulate random realisations (random intervals) from the inverses of the bounding CDFs

$$     \Gamma(\alpha)=[\overline{F}^{-1}(\alpha), \underline{F}^{-1}(\alpha)],$$

where $\alpha \sim U(0,1)$ is a sample from a uniform distribution. This interval can be visualised as a horizontal cut (or slice) of the p-box at an $\alpha \in [0,1]$. This interval is then propagated through the model $f$ or performance function $g$ using an optimiser or surrogate model,

$$     \underline{g}(\alpha) = \min_{x \in \Gamma(\alpha)}g(x),$$

$$
     \overline{g}(\alpha) = \max_{x \in \Gamma(\alpha)}g(x).$$

In _UncertaintyQuantification.jl_ the intervals are also propagated using MADS. In the multivariate case, we can combine two correlated intervals using a Cartesian product

$$     [\overline{F}_X^{-1}(\alpha_X), \underline{F}_X^{-1}(\alpha_X)] \times [\overline{F}_Y^{-1}(\alpha_Y), \underline{F}_Y^{-1}(\alpha_Y)],$$

where $(\alpha_X, \alpha_Y) \sim C_{XY}$ are samples of the copula between $X$ and $Y$.

The reliability analysis can be written as thus:

$$     \underline{p}*{\text{f}} = \int_U \mathbb{I}[\overline{g}(\alpha)] dC,(\alpha)$$

$$     \overline{p}*{\text{f}} = \int_U \mathbb{I}[\underline{g}(\alpha)] dC.(\alpha)$$

In some sense, the two loops from the double-loop method have been reversed, where now the outer-loop handles the random (aleatory) component, and the inner-loop handles the interval propagation (epistemic). Describing the analysis this way essentially gives two separate reliability calculations, with $\underline{g}$ and $\overline{g}$ as the two target performance functions. Rosenblatt transformations may be used to associate a standard normal distribution to the copula $C$, and one may then use any standard reliability method to compute the performance bounds.

The software implementation is such that this imprecise reliability method can be coupled to any simulation method (FORM, line sampling, etc.) in a straightforward way. As with the double-loop, one can apply different simulations for each bound if desired.

The problem setup for random slicing is identical to that of the double-loop. The only difference is that a [`RandomSlicing`](/api/reliability#UncertaintyQuantification.RandomSlicing) instance is passed instead of [`DoubleLoop](/api/reliability#UncertaintyQuantification.DoubleLoop).

Here, the lower bound is again estimated using FORM, while we apply subset simulation to obtain the upper bound. The result for the lower bound matches the analytical solution perfectly. The upper bound is estimated accurately as $\overline{p}_f \approx 3.884325e-5$. Note, that in addition to the probability of failure, random slicing also returns other outputs of the underlying simulations, such as the coefficient of variation and the evaluated samples for potential post-processing.

```julia
     subset = SubSetSimulation(2000, 0.1, 10, Normal())
     pf, out_lb, out_ub = probability_of_failure(f, df -> 9 .+ df.f, [x1,x2], RandomSlicing(FORM(), subset))
```


```ansi
([7.619853024348675e-24, 2.0704482375000004e-5], (9.999999999997554, (x2 = -9.999999999997554,), (x2 = 1.0,)), (5.681608394188276e-6, 10000×3 DataFrame
   Row │ x2                                 g_slice     level
       │ Interval                           Float64     Int64
───────┼──────────────────────────────────────────────────────
     1 │ [-1.856648965655265, -0.92832448…   6.14335        1
     2 │ [-0.30592702592635285, -0.152963…   7.69407        1
     3 │ [-0.10477136303960666, -0.052385…   7.89523        1
     4 │ [-2.6971631917656045, -1.3485815…   5.30284        1
     5 │ [-0.2326571514017165, -0.1163285…   7.76734        1
     6 │ [-0.4206885399014404, -0.2103442…   7.57931        1
     7 │ [-0.0981693778560777, -0.0490846…   7.90183        1
     8 │ [0.3534422809414995, 0.706884561…   8.35344        1
   ⋮   │                 ⋮                      ⋮         ⋮
  9994 │ [-7.360701731357478, -3.68035086…   0.639298       5
  9995 │ [-7.2516490572991925, -3.6258245…   0.748351       5
  9996 │ [-7.246100087534909, -3.62305004…   0.7539         5
  9997 │ [-8.739034238143178, -4.36951711…  -0.739034       5
  9998 │ [-7.386167369040863, -3.69308368…   0.613833       5
  9999 │ [-8.50239659615072, -4.251198298…  -0.502397       5
 10000 │ [-8.059425665032638, -4.02971283…  -0.0594257      5
                                             9985 rows omitted))
```


As outlined by [alvarez2018estimation], other forms of random-sets can in principle be evaluated with this method, such as possibility distributions [[22](/references#dubois1990consonant)] or general Dempster–Shafer structures [[23](/references#shafer1976mathematical)]. However, careful consideration of multivariate extensions of these structures must be taken [[24](/references#schmelzer2023random)]. For this reason, we restrict ourselves to distributions, intervals, and p-boxes for the time being.
