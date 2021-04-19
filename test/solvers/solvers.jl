@testset "Solvers" begin

    binary = ""
    if Sys.iswindows()
    binary = joinpath(pwd(), "solvers/bin/radius.exe")
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
end
