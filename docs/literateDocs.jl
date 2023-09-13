using Literate

function fuseConvert(path::String, out::IOStream, r::String, dir::String)
    for (root, _, files) in walkdir(joinpath(r, dir))
        for file in files
            write(path, read(joinpath(root, file), String) * "\n\n")
        end
    end
    close(out)
    Literate.markdown(path, "./docs/src/examples"; documenter=true, name=dir)

    #remove @meta block created by literate
    s = read("./docs/src/examples/" * dir * ".md", String)
    lines = split(s, "\n"; limit=4)
    open("./docs/src/examples/" * dir * ".md", "w") do f
        write(f, lines[4])
        close(f)
    end

    return nothing
end

for (r, d, f) in walkdir("./docs/literate/")
    for dir in d
        tmpPath, tmpFile = mktemp()
        fuseConvert(tmpPath, tmpFile, r, dir)
    end
end

# out = ""
#         for (root, dirs, files) in walkdir(joinpath(r, dir))
#             for file in files
#                 if (cmp(file, dir * ".jl") != 0)
#                     out = out * read(joinpath(root, file), String) * "\n\n"
#                 end
#             end
#         end

#         open(joinpath(joinpath(r, dir), dir * ".jl"), "w") do f
#             write(f, out)
#             close(f)
#         end

#         Literate.markdown(
#             joinpath(joinpath(r, dir), dir * ".jl"), "./docs/src/examples"; documenter=true
#         )
