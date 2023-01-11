using UncertaintyQuantification

x = RandomVariable.(Normal(), [:x1, :x2, :x3, :x4])

function limitstate(df)
    return exp.(-0.3671 * (df.x1 .+ 2df.x2 .+ 3df.x3)) .- df.x4 .+ 1.5
end

@time pf, β, dp = probability_of_failure(limitstate, x, FORM())
