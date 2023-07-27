using Literate

for (root, dirs, files) in walkdir("./docs/literate/")
    if (!isempty(files))
        Literate.markdown.(joinpath.(root, files), "./docs/src/examples", documenter=true)

        sp = splitpath(root)
        if (cmp(sp[end], "literate") == 0)
            Literate.script.(joinpath.(root, files), joinpath.(pwd(), "./demo"))
        else
            Literate.script.(
                joinpath.(root, files), joinpath.(pwd(), "./demo/", joinpath(sp[4:end]))
            )
        end
    end
end
