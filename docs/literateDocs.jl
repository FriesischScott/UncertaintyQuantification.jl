using Literate

"""
    function fuseConvert(path::String, out::IOStream, r::String, dir::String)

Groups all files from a directory into a temp file and creates one markdown-file (saved to ./docs/src/examples) using Literate.jl.

"""
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