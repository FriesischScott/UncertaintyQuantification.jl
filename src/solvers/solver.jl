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

    out = joinpath(folder, basename(binary) * "UncertaintyQuantification.out")
    err = joinpath(folder, basename(binary) * "UncertaintyQuantification.err")

    p = pipeline(
        !isempty(args) ? `$binary $args $source` : `$binary $source`; stdout=out, stderr=err
    )

    cd(() -> run(p), folder)

    return nothing
end
