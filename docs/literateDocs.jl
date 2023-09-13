using Literate

function fuseConvert(path::String, out::IOStream, r::String, dir::String)
    for (root, dirs, files) in walkdir(joinpath(r, dir))
        for file in files
            write(path, read(joinpath(root, file), String) * "\n\n")
        end
    end
    close(out)
    Literate.markdown(path, "./docs/src/examples"; documenter=true, name=dir)
    return nothing
end

for (r, d, f) in walkdir("./docs/literate/")
    for dir in d
        tmpPath, tmpFile = mktemp()
        fuseConvert(tmpPath, tmpFile, r, dir)
    end
end
