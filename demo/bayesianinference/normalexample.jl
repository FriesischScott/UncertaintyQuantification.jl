using UncertaintyQuantification
using DataFrames
using Plots

nSamples = 2000
nChains = 16

x0 = DataFrame([:x, :y, :z] .=> [0.0, 0.0, 0.0])
C0 = [[1.,0.,0.] [0.,1.,0.] [0.,0.,1.]]

amhSampler = AdaptiveMetropolisHastings(nSamples,nChains,x0,C0)

prior = Normal();
priorF = df -> pdf.(prior,df.x) .* pdf.(prior,df.y) .* pdf.(prior,df.z)

logL = df -> -((df.x.+3).^2 .+ (df.y.-2).^2 .+ (df.z.+1).^2)./1e-2;

samples = bayesianupdating(logL,priorF,amhSampler);

print(size(samples))

scatter3d(samples.x,samples.y,samples.z)
