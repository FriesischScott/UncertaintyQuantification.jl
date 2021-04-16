struct Solver
    path::String
    args::String
    source::String
end

function run(solver::Solver, folder::String)
    binary = solver.path
    args = solver.args
    source = solver.source

    old_pwd = pwd()
    cd(folder)

    out = joinpath(folder, basename(binary) * "UncertaintyQuantification.out")
    err = joinpath(folder, basename(binary) * "UncertaintyQuantification.err")

    if !isempty(args)
        Base.run(pipeline(`$binary $args $source`, stdout=out, stderr=err))
    else
        Base.run(pipeline(`$binary $source`, stdout=out, stderr=err))
    end

    cd(old_pwd)
end