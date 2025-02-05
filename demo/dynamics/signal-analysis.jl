#using UncertaintyQuantification

# Parallel processing of MC simulation
using Distributed
addprocs(6; exeflags="--project")
@everywhere using UncertaintyQuantification 

ω = collect(0:0.6:150)

cp = CloughPenzien(ω, 0.1, 0.8π, 0.6, 8π, 0.6)

# Ground motion model
gm = SpectralRepresentation(cp, collect(0:0.02:10), :gm)
gm_model = StochasticProcessModel(gm)

capacity = Parameter(21, :cap) #estimation for a capacity value resulting in pf ~ [1.0e-5, 1.0e-6]

function limitstate(df)
    return df.cap - map(x -> maximum(abs.(x)), df.gm)
end

models = [gm_model]
inputs = [gm, capacity]

# Compute probability of failure using standard Monte Carlo
mc = MonteCarlo(10^4)

@time mc_pf, mc_std, mc_samples = probability_of_failure(models, limitstate, inputs, mc)

println(
    "Monte Carlo probability of failure $mc_pf ($(size(mc_samples, 1)) model evaluations)"
)

# Compute probability of failure using Importance Sampling
pf, β, dp, α = probability_of_failure(models, limitstate, inputs, FORM())
is = ImportanceSampling(10^4, β, dp, α)

is_pf, is_std, is_samples = probability_of_failure(models, limitstate, inputs, is)

println(
    "Importance Sampling probability of failure: $is_pf ($(size(is_samples, 1)) model evaluations)",
)

# Compute probability of failure using Line Sampling
ls = LineSampling(200)

ls_pf, ls_std, ls_samples = probability_of_failure(models, limitstate, inputs, ls)

println(
    "Line Sampling probability of failure: $ls_pf ($(size(ls_samples, 1)) model evaluations)",
)

# Compute probability of failure using Subset Sampling
subset = UncertaintyQuantification.SubSetSimulation(2000, 0.1, 10, Uniform(-0.5, 0.5))

subset_pf, subset_std, subset_samples = probability_of_failure(models, limitstate, inputs, subset)

println(
    "Subset Simulation probability of failure: $subset_pf ($(size(subset_samples, 1)) model evaluations)",
)

# Compute probability of failure using conditional Subset Sampling
subset_inf = UncertaintyQuantification.SubSetInfinity(2000, 0.1, 10, 0.5)

subset_pf_inf, subset_std_inf, subset_samples_inf = probability_of_failure(models, limitstate, inputs, subset_inf)

println(
    "Subset infinity probability of failure: $subset_pf_inf ($(size(subset_samples_inf, 1)) model evaluations)",
)

# Compute probability of failure using adaptive Subset Sampling
subset_adap = UncertaintyQuantification.SubSetInfinityAdaptive(2000, 0.1, 10, 10, 0.6, 1.0)

subset_pf_adap, subset_std_adap, subset_samples_adap = probability_of_failure(models, limitstate, inputs, subset_adap)

println(
    "Subset infinity adaptive probability of failure: $subset_pf_adap ($(size(subset_samples_adap, 1)) model evaluations)",
)