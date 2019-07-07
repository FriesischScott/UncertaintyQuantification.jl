module OpenCossan

using Distributions, LinearAlgebra, DataFrames

import Base: rand

export
    # inputs
    RandomVariableSet

    # methods
    rand

include("inputs/randomvariableset.jl")

end