struct Solver
    path::String
    source::String
    args::String
end

Solver(path::String, source::String; args::String="") = Solver(path, source, args)

function run(solver::Solver, folder::String)
    binary = solver.path
    source = solver.source
    args = solver.args

    old_pwd = pwd()
    cd(folder)

    out = joinpath(folder, basename(binary) * "UncertaintyQuantification.out")
    err = joinpath(folder, basename(binary) * "UncertaintyQuantification.err")

    if !isempty(args)
        Base.run(pipeline(`$binary $args $source`; stdout=out, stderr=err))
    else
        Base.run(pipeline(`$binary $source`; stdout=out, stderr=err))
    end

    return cd(old_pwd)
end
