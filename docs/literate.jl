using Literate

for (root, dirs, files) in walkdir("./docs/literate/")
    if (!isempty(files))
        Literate.markdown.(joinpath.(root, files), "./docs/src/examples", documenter=true)

        sp = splitpath(root)

        Literate.script.(
            joinpath.(root, files),
            joinpath.(
                pwd(),
                "./demo/",
                joinpath(sp[(findfirst(item -> item == "literate", sp) + 1):end]),
            ),
        )
    end
end
