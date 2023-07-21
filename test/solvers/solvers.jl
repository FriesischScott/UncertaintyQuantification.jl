@testset "Solvers" begin
    binary = ""
    if Sys.iswindows()
        binary = joinpath(pwd(), "solvers/bin/radius.exe")
    elseif Sys.isapple()
        binary = joinpath(pwd(), "solvers/bin/radius-mac")
    else
        binary = joinpath(pwd(), "solvers/bin/radius")
    end

    solver = Solver(binary, "", "in.txt")

    tmp = tempdir()

    open(joinpath(tmp, "in.txt"), "w") do input
        println(input, "0.5")
        println(input, "0.5")
    end

    run(solver, tmp)

    radius = map(x -> parse(Float64, x), readlines(joinpath(tmp, "out.txt")))

    @test radius[1] == sqrt(0.5^2 + 0.5^2)

    current_dir = pwd()

    @test_logs (
        :warn,
        "Solver path is not absolute. Make sure this_path_is_not_absolute is on your PATH.",
    ) Solver("this_path_is_not_absolute", "in.txt")

    @test_throws Base.IOError run(
        Solver(joinpath(pwd(), "this/solver/does/not/exist"), "in.txt"), tmp
    )
    @test pwd() == current_dir
end
