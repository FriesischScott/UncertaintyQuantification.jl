using Literate

for (root, dirs, files) in walkdir("./docs/literate/")
    if (!isempty(files))
        Literate.markdown.(joinpath.(root, files), "./docs/src/examples", documenter=true)
    end
end
