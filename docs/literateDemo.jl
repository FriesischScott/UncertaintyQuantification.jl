using Literate

for (root, dirs, files) in walkdir("./docs/literate/")
    sp = splitpath(root)
    if (!isempty(files))
        Literate.script.(
            joinpath.(root, files), joinpath.(pwd(), "./demo/", joinpath(sp[4:end]))
        )
    end
end