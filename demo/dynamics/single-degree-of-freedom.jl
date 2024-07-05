using UncertaintyQuantification
using Dierckx

using DifferentialEquations

m = Parameter(0.5, :m)

k = Parameter(2, :k)
c = Parameter(0, :c)

ω = collect(0:0.6:150)

cp = CloughPenzien(ω, 0.1, 0.8π, 0.6, 8π, 0.6)

ex = SpectralRepresentation(cp, collect(0:0.02:10), :ex)

function sdof(df)
    return map(eachrow(df)) do s
        excitation = Spline1D(ex.time, ex(s); k=1)

        function f(dy, y, _, t)
            dy[1] = y[2]

            return dy[2] = -s.c / s.m * y[2] - s.k / s.m * y[1] + excitation(t)
        end

        prob = ODEProblem(f, [0.0, 0.0], (ex.time[1], ex.time[end]))

        sol = solve(prob, Tsit5())

        return sol[1, :]
    end
end

displacement = Model(sdof, :d)

inputs = [ex, m, k, c]

# samples = sample(inputs, 10)

# evaluate!(displacement, samples)

# p = plot()

# for i in 1:size(samples, 1)
#     t = range(ex.time[1], ex.time[end], length(samples[i, :d]))

#     plot!(p, t, samples[i, :d])
# end

function g(df)
    return map(eachrow(df)) do s
        1.5 - maximum(abs.(s.d))
    end
end

mc_pf, mc_std, mc_samples, = probability_of_failure(
    [displacement], g, inputs, MonteCarlo(10000)
)

ss_pf, ss_std, ss_samples = probability_of_failure(
    [displacement], g, inputs, SubSetSimulation(2000, 0.1, 10, Uniform(-0.5, 0.5))
)
