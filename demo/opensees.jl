using UncertaintyQuantification, DelimitedFiles, Formatting

dt = Parameter(0.001, :dt)
steps = Parameter(35 / dt.value, :steps)

s1p = RandomVariable(Normal(10, 5), :s1p)

uni1 = RandomVariable(Uniform(5, 10), :uni1)
uni2 = RandomVariable(Uniform(0, 0.5), :uni2)

e1p = RandomVariable(Uniform(0.0002, 0.0003), :e1p)
e2p = RandomVariable(Uniform(0.0006, 0.0009), :e2p)
e3p = RandomVariable(Uniform(0.0013, 0.0020), :e3p)

e1n = Model(df -> -1 .* df.e1p, :e1n)
e2n = Model(df -> -1 .* df.e2p, :e2n)
e3n = Model(df -> -1 .* df.e3p, :e3n)

s2p = Model(df -> df.s1p .+ df.uni1, :s2p)
s3p = Model(df -> df.s1p, :s3p)

s1n = Model(df -> -1 .* df.s1p, :s1n)
s2n = Model(df -> -1 .* df.s2p, :s2n)
s3n = Model(df -> df.s2n .+ df.s2n .* -1 .* df.uni2, :s3n)

inputs = [dt, steps, s1p, uni1, uni2, e1p, e2p, e3p]

# Source/Extra files are expected to be in this folder
sourcedir = "/home/behrensd/Seafile/Documents/mz-opensees"

# These files will be rendere through Mustach.jl and have values injected
sourcefiles = [
    "0.SteelFrame.tcl",
    "5.Materials.tcl"
]

# These files will be copied to the working directory without injecting values
extrafiles = [
    "1.NodeCoord.tcl",
    "2.SPConstraint.tcl",
    "3.NodeMass.tcl",
    "4.MPConstraint.tcl",
    "6.Sections.tcl",
    "7.GeoTran.tcl",
    "8.Elements.tcl",
    "9.TimeSeries.tcl",
    "10.LoadPattern_1.tcl",
    "11.Recorder_1.tcl",
    "12.AnalysisOptn_1.tcl",
    "13.LoadPattern_2.tcl",
    "14.Recorder_2.tcl",
    "15.AnalysisOptn_2.tcl",
    "LibAnalysisDynamicParameters.tcl",
    "record-kozani-mainshock.txt"
]

# Dictionary to map FormatSpecs (Formatting.jl) to variables
numberformats = Dict(
    :dt => FormatSpec("f"),
    :steps => FormatSpec("d"),
    :* => FormatSpec("16.8e") # Catches all other variables
)

# UQ will create subfolders in here to run the solver and store the results
workdir = "/home/behrensd/src/github.com/friesischscott/UncertaintyQuantification.jl/opensees"

# Read output file and compute maximum (absolute) displacement
# An extractor is based the working directory for the current sample
d_max = Extractor(
    base -> begin
        file = joinpath(base, "ef_kozms_Dsp.txt")
        data = readdlm(file, ' ')

        return maximum(abs.(data[:, 2]))
    end,
    :d_max
)

# Multiple extractors can be passed to the model
# Read output file and compute minimum (absolute) displacement
d_min = Extractor(
    base -> begin
        file = joinpath(base, "ef_kozms_Dsp.txt")
        data = readdlm(file, ' ')

        return minimum(abs.(data[:, 2]))
    end,
    :d_min
)

extractors = [d_max, d_min]

opensees = Solver(
    "/home/behrensd/src/OpenSees/bin/OpenSees",
    "",
    "0.SteelFrame.tcl"
)

ext = ExternalModel(
    sourcedir,
    sourcefiles,
    extrafiles,
    numberformats,
    workdir,
    extractors,
    opensees
)

models = [e1n, e2n, e3n, s2p, s3p, s1n, s2n, s3n, ext]

df = sample(inputs, 2)

evaluate!(models, df)

println("Maximum displacement: $(df.d_max)")
println("Minimum displacement: $(df.d_min)")