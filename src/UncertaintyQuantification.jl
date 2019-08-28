module UncertaintyQuantification

using Distributions, LinearAlgebra, DataFrames

import Base: rand

export
    # inputs
    RandomVariableSet,

    # methods
    copularand,
    rand

include("inputs/randomvariableset.jl")
include("util/copula.jl")

end