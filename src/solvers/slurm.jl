## Probably better to make a 'SlurmSolver' version of 'solver'

struct SlurmInterface

    name::String
    account::String
    nodes::Integer
    ntasks::Integer
    time::String
    partition::String
    extras::Vector{String}

end

function make_input(SlurmInterface)

end

function launch(SlurmInterface, ExternalModel )
    
end