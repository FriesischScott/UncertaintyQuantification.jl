
module TypstWriter
using Documenter
using MarkdownAST: MarkdownAST, Node
using AbstractTrees

struct Typst <: Documenter.Writer
    platform::String
    version::String
    bib::String
    authors::Vector{String}
    affiliations::Vector{String}
    function TypstWriter.Typst(;
        platform="native",
        version=get(ENV, "TRAVIS_TAG", ""),
        bib="",
        authors=String[],
        affiliations=String[],
    )
        platform ∈ ("native", "docker", "none") ||
            throw(ArgumentError("unknown platform: $platform"))
        return new(platform, string(version), bib, authors, affiliations)
    end
end

abstract type TypstFormat <: Documenter.FormatSelector end
Documenter.Selectors.order(::Type{TypstFormat}) = 4.0
Documenter.Selectors.matcher(::Type{TypstFormat}, fmt, _) = isa(fmt, TypstWriter.Typst)
Documenter.Selectors.runner(::Type{TypstFormat}, fmt, doc) = TypstWriter.render(doc, fmt)

abstract type CitationPipeline <: Documenter.Builder.DocumentPipeline end

Documenter.Selectors.order(::Type{CitationPipeline}) = 2.1

struct Citation <: MarkdownAST.AbstractElement
    form::String
    key::String
end

const cite_map = Dict{String,String}("@cite" => "normal", "@citet" => "prose")

function Citation(node::Node)
    if !haskey(cite_map, node.element.destination)
        error("Unknown citation type $(node.element.destination)")
    end
    return Citation(cite_map[node.element.destination], first(node.children).element.text)
end

function Documenter.Selectors.runner(_::Type{CitationPipeline}, doc::Documenter.Document)
    @info "Processing citations."
    for (_, page) in doc.blueprint.pages
        for node in AbstractTrees.PreOrderDFS(page.mdast)
            if node.element isa MarkdownAST.Link &&
                occursin("@cite", node.element.destination)
                node.element = Citation(node)
                empty!(node.children)
            end
            if node.element isa Documenter.DocsNode
                for mdast in node.element.mdasts
                    for n in AbstractTrees.PreOrderDFS(mdast)
                        if n.element isa MarkdownAST.Link &&
                            occursin("@cite", n.element.destination)
                            n.element = Citation(n)
                            empty!(n.children)
                        end
                    end
                end
            end
        end
    end
    return nothing
end

using Documenter: Documenter
using Markdown: Markdown
using ANSIColoredPrinters: ANSIColoredPrinters

mutable struct Context{I<:IO} <: IO
    io::I
    in_header::Bool
    footnotes::Dict{String,Int}
    depth::Int
    filename::String # currently active source file
    bib::String
    doc::Documenter.Document
end

Context(io, doc) = Context{typeof(io)}(io, false, Dict(), 1, "", "", doc)

_print(c::Context, args...) = Base.print(c.io, args...)
_println(c::Context, args...) = Base.println(c.io, args...)
_print(io, args...) = Base.print(io, args...)

# Labels in the TeX file are hashes of plain text labels.
# To keep the plain text label (for debugging), say _hash(x) = x
_hash(x) = string(hash(x))

const STYLE = joinpath(dirname(@__FILE__), "..", "..", "assets", "latex", "documenter.sty")
const DEFAULT_PREAMBLE_PATH = joinpath(
    dirname(@__FILE__), "..", "..", "assets", "latex", "preamble.tex"
)

function hastypst()
    try
        return success(`typst --version`)
    catch
        return false
    end
end

const DOCUMENT_STRUCTURE = ("=", "==", "===", "====", "=====")

function render(doc::Documenter.Document, settings::Typst=Typst())
    if isempty(doc.user.sitename) # otherwise, the latex compiler will terminate with a cryptic "There's no line here to end" error
        error(
            """
            LaTeXWriter needs a non-empty `sitename` passed to `makedocs`, otherwise the LaTeX build will error!
            Please pass e.g. `sitename = "Some Site Name"` as a keyword argument to `makedocs`.
            """,
        )
    end

    @info "TypstWriter: creating the Typst file."
    mktempdir() do path
        cp(joinpath(doc.user.root, doc.user.build), joinpath(path, "build"))
        if !isempty(settings.bib)
            cp(settings.bib, joinpath(path, "build", basename(settings.bib)))
        end
        cd(joinpath(path, "build")) do
            fileprefix = latex_fileprefix(doc, settings)
            open("$(fileprefix).typ", "w") do io
                context = Context(io, doc)
                context.bib = basename(settings.bib)
                writeheader(context, doc, settings)
                for (title, filename, depth) in files(doc.user.pages)
                    context.filename = filename
                    context.depth = depth - 1
                    empty!(context.footnotes)
                    if isempty(filename) && depth == 1
                        _println(context, "$(DOCUMENT_STRUCTURE[1]) $title")
                    else
                        path = normpath(filename)
                        page = doc.blueprint.pages[path]
                        if get(page.globals.meta, :IgnorePage, :none) !== :typst
                            typst(context, page.mdast.children; toplevel=true)
                        end
                    end
                end
            end
            # compile .tex
            status = compile_tex(doc, settings, fileprefix)
            # Debug: if DOCUMENTER_LATEX_DEBUG environment variable is set, copy the LaTeX
            # source files over to a directory under doc.user.root.
            if haskey(ENV, "DOCUMENTER_TYPST_DEBUG")
                dst = if isempty(ENV["DOCUMENTER_TYPST_DEBUG"])
                    mktempdir(doc.user.root; cleanup=false)
                else
                    joinpath(doc.user.root, ENV["DOCUMENTER_TYPST_DEBUG"])
                end
                sources = cp(pwd(), dst; force=true)
                @info "Typst sources copied for debugging to $(sources)"
            end

            # If the build was successful, copy the PDF or the LaTeX source to the .build directory
            if status && (settings.platform != "none")
                pdffile = "$(fileprefix).pdf"
                cp(pdffile, joinpath(doc.user.root, doc.user.build, pdffile); force=true)
            elseif status && (settings.platform == "none")
                cp(pwd(), joinpath(doc.user.root, doc.user.build); force=true)
            else
                error("Compiling the .typ file failed. See logs for more information.")
            end
        end
    end
    return nothing
end

function latex_fileprefix(doc::Documenter.Document, settings::Typst)
    fileprefix = doc.user.sitename
    if occursin(Base.VERSION_REGEX, settings.version)
        v = VersionNumber(settings.version)
        fileprefix *= "-$(v.major).$(v.minor).$(v.patch)"
    end
    return replace(fileprefix, " " => "")
end

const DOCKER_IMAGE_TAG = "0.1"

function compile_tex(doc::Documenter.Document, settings::Typst, fileprefix::String)
    if settings.platform == "native"
        Sys.which("typst") === nothing &&
            (@error "TypstWriter: typst command not found."; return false)
        @info "TypstWriter: using typst to compile."
        try
            piperun(`typst compile $(fileprefix).typ`; clearlogs=true)
            return true
        catch err
            logs = cp(pwd(), mktempdir(; cleanup=false); force=true)
            @error "LaTeXWriter: failed to compile tex with latexmk. " *
                "Logs and partial output can be found in $(Documenter.locrepr(logs))" exception =
                err
            return false
        end
    elseif settings.platform == "docker"
        Sys.which("docker") === nothing &&
            (@error "LaTeXWriter: docker command not found."; return false)
        @info "LaTeXWriter: using docker to compile tex."
        script = """
        mkdir /home/zeptodoctor/build
        cd /home/zeptodoctor/build
        cp -r /mnt/. .
        latexmk -f -interaction=batchmode -halt-on-error -view=none -lualatex -shell-escape $(fileprefix).tex
        """
        try
            piperun(
                `docker run -itd -u zeptodoctor --name latex-container -v $(pwd()):/mnt/ --rm juliadocs/documenter-latex:$(DOCKER_IMAGE_TAG)`;
                clearlogs=true,
            )
            piperun(`docker exec -u zeptodoctor latex-container bash -c $(script)`)
            piperun(`docker cp latex-container:/home/zeptodoctor/build/$(fileprefix).pdf .`)
            return true
        catch err
            logs = cp(pwd(), mktempdir(; cleanup=false); force=true)
            @error "LaTeXWriter: failed to compile tex with docker. " *
                "Logs and partial output can be found in $(Documenter.locrepr(logs))" exception =
                err
            return false
        finally
            try
                piperun(`docker stop latex-container`)
            catch
            end
        end
    elseif settings.platform == "none"
        @info "Skipping compiling typst file."
        return true
    end
end

function piperun(cmd; clearlogs=false)
    verbose = "--verbose" in ARGS || get(ENV, "DOCUMENTER_VERBOSE", "false") == "true"
    if verbose
        cmd = pipeline(
            cmd; stdout="LaTeXWriter.stdout", stderr="LaTeXWriter.stderr", append=!clearlogs
        )
    end
    return run(cmd)
end

function writeheader(io::IO, doc::Documenter.Document, settings::Typst)
    # custom = joinpath(doc.user.root, doc.user.source, "assets", "custom.sty")
    # isfile(custom) ? cp(custom, "custom.sty"; force=true) : touch("custom.sty")

    # custom_preamble_file = joinpath(doc.user.root, doc.user.source, "assets", "preamble.tex")
    # if isfile(custom_preamble_file)
    #     # copy custom preamble.
    #     cp(custom_preamble_file, "preamble.tex"; force=true)
    # else # no custom preamble.tex, use default.
    #     cp(DEFAULT_PREAMBLE_PATH, "preamble.tex"; force=true)
    # end
    preamble = """
#import "@preview/mitex:0.2.5": *
#import "@preview/note-me:0.5.0": *

#set page(
  paper: "a4"
)

#show heading: set block(below: 2em, above: 2em)
#show heading: it => align(left)[#it]

#show raw.where(block:true): it => block(
  fill: rgb("#8e96aa24"),
  inset: 8pt,
  radius: 4pt,
  text(font:"Fira Code", it)
)

#show raw.where(block:false): it => box(
  fill: rgb("#8e96aa24"),
  inset: (x: 3pt, y: 0pt),
  outset: (y: 3pt),
  radius: 2pt,
  text(font:"Fira Code", fill: rgb("#3451b2"), it)
)

#let title = "$(doc.user.sitename)"
#let version = "$(doc.user.version)"

#set page(numbering: {})

#set align(center)
#text(weight: "bold", size: 20pt, title)
#v(0.5em)
#text(style: "italic", size: 12pt, version)

#par(justify: true)[
    *Abstract* \
    #lorem(50)
  ]

#set page(numbering: "I")

#set par(justify: true)
#set text(
  font: "New Computer Modern",
  size: 12pt,
)

#outline()

#counter(page).update(1)
#set page(numbering: "1")

#set heading(numbering: "1.1.1")

"""
    _println(io, preamble)
    return nothing
end

# A few of the nodes are printed differently depending on whether they appear
# as the top-level blocks of a page, or somewhere deeper in the AST.
istoplevel(n::Node) = !isnothing(n.parent) && isa(n.parent.element, MarkdownAST.Document)

typst(io::Context, node::Node) = typst(io, node, node.element)

function typst(io::Context, node::Node, e::MarkdownAST.AbstractElement)
    @warn "$(typeof(e)) not implemented: $e"
end

function typst(io::Context, children; toplevel=false)
    @assert eltype(children) <: MarkdownAST.Node
    for node in children
        # otherelement = !isa(node.element, NoExtraTopLevelNewlines)
        # toplevel && otherelement && _println(io)
        typst(io, node)
        # toplevel && otherelement && _println(io)
    end
    return nothing
end

function typst(io::Context, node::Node, mcb::Documenter.MultiCodeBlock)
    @show mcb
    return nothing
end

using Base64: base64decode

function typst(io::Context, node::Node, heading::MarkdownAST.Heading)
    _print(io, "$(DOCUMENT_STRUCTURE[heading.level+io.depth]) ")
    typst(io, node.children)
    _println(io, "")
    return nothing
end

function typst(io::Context, node::Node, list::MarkdownAST.List)
    @show list
    return nothing
end

function typst(io::Context, node::Node, code::MarkdownAST.CodeBlock)
    language = Documenter.codelang(code.info)
    if language == "@bibliography"
        _println(io, "#bibliography(title: none, \"$(io.bib)\")")
    else
        text = IOBuffer(code.code)
        code_code = repr(MIME"text/plain"(), ANSIColoredPrinters.PlainTextPrinter(text))
        _println(io, "```$(language)")
        _println(io, code_code)
        _println(io, "```\n")
    end
    return nothing
end

function typst(io::Context, _::Node, text::MarkdownAST.Text)
    latexesc(io, text.text)
    return nothing
end

function typst(io::Context, node::Node, code::MarkdownAST.Code)
    _print(io, "`$(code.code)`")
    return nothing
end

function typst(io::Context, node::Node, _::MarkdownAST.Strong)
    _print(io, "*")
    typst(io, node.children)
    _print(io, "*")
    return nothing
end

function typst(io::Context, node::Node, _::MarkdownAST.Emph)
    _print(io, "_")
    typst(io, node.children)
    _print(io, "_")
    return nothing
end

function typst(io::Context, node::Node, link::MarkdownAST.Link)
    if !isempty(node.children)
        _print(io, "#link(\"$(link.destination)\")[")
        typst(io, node.children)
        _print(io, "]")
    else
        _print(io, "#link(\"$(link.destination)\")")
    end
    return nothing
end

function typst(io::Context, node::Node, math::MarkdownAST.InlineMath)
    _print(io, "#mi(\"$(math.math)\")")
    return nothing
end

function typst(io::Context, node::Node, math::MarkdownAST.DisplayMath)
    _println(io, "#mitex(`")
    _print(io, math.math)
    _println(io, "`)\n")
    return nothing
end

function typst(io::Context, node::Node, head::Documenter.AnchoredHeader)
    typst(io, node.children)
    return nothing
end

function typst(io::Context, ::Node, d::Dict{MIME,Any})
    filename = String(rand('a':'z', 7))
    if haskey(d, MIME"image/png"())
        @show "png"
        write("$(filename).png", base64decode(d[MIME"image/png"()]))
        # _println(
        #     io, """
        #     \\begin{figure}[H]
        #     \\centering
        #     \\includegraphics[max width=\\linewidth]{$(filename)}
        #     \\end{figure}
        #     """
        # )
    elseif haskey(d, MIME"image/jpeg"())
        write("$(filename).jpeg", base64decode(d[MIME"image/jpeg"()]))
        # _println(
        #     io, """
        #     \\begin{figure}[H]
        #     \\centering
        #     \\includegraphics[max width=\\linewidth]{$(filename)}
        #     \\end{figure}
        #     """
        # )
        @show "png"
    elseif haskey(d, MIME"text/latex"())
        @show "latex"
        # If it has a latex MIME, just write it out directly.
        # content = d[MIME("text/latex")]
        # if startswith(content, "\\begin{tabular}")
        #     # This is a hacky fix for the printing of DataFrames (or any type
        #     # that produces a {tabular} environment). The biggest problem is
        #     # that tables with may columns will run off the page. An ideal fix
        #     # would be for the printing to omit some columns, but we don't have
        #     # the luxury here. So instead we just rescale everything until it
        #     # fits. This might make the rows too small, but it's arguably better
        #     # than having them go off the page.
        #     _println(io, "\\begin{table}[h]\n\\centering")
        #     _println(io, "\\adjustbox{max width=\\linewidth}{")
        #     _print(io, content)
        #     _println(io, "}\\end{table}")
        # else
        # _print(io, content)
        # end
    elseif haskey(d, MIME"text/markdown"())
        @show "md"
        # md = Markdown.parse(d[MIME"text/markdown"()])
        # ast = MarkdownAST.convert(MarkdownAST.Node, md)
        # latex(io, ast.children)
    elseif haskey(d, MIME"text/plain"())
        text = d[MIME"text/plain"()]
        out = repr(MIME"text/plain"(), ANSIColoredPrinters.PlainTextPrinter(IOBuffer(text)))

        codeblock = MarkdownAST.CodeBlock("julia", out)
        typst(io, MarkdownAST.Node(codeblock))
    else
        error("this should never happen.")
    end
    return nothing
end

function typst(io::Context, node::Node, multi::Documenter.MultiOutput)
    return typst(io, node.children)
end

function typst(io::Context, node::Node, multi::Documenter.MultiOutputElement)
    return typst(io, node, multi.element)
end

function typst(io::Context, node::Node, img::Documenter.LocalImage)
    _print(io, "#image(\"$(replace(img.path,"\\" => "/"))\")")
    return nothing
end

function typst(io::Context, node::Node, cite::Citation)
    keys = split(cite.key, ',')
    for (i, k) in enumerate(keys)
        _print(io, "#cite(<$(k)>, form: \"$(cite.form)\")")
        if i < length(keys)
            _print(io, " ")
        end
    end
end

function typst(io::Context, node::Node, docs::Documenter.DocsNodesBlock)
    typst(io, node.children)
    return nothing
end

function typst(io::Context, node::Node, docs::Documenter.DocsNode)
    # TODO: Add link
    typst(io, node.children)
    return nothing
end

function typst(io::Context, node::Node, index::Documenter.IndexNode)
    for el in index.elements
        @show el
    end
end

# function typst(io::Context, node::Node, link::Documenter.PageLink)
#     # if length(node.children) >0 && isa(node.children[1], MarkdownAST.Text) && isa(tryparse(Int, node.children[1].text), Integer)
#     #     @show "Found Citation pointing to $(link.fragment)"
#     # end
#     @show node
#     @show is_citation(node)
# end

# typst(io::Context, node::Node, mcb::Documenter.MultiCodeBlock) = typst(io, node, join_multiblock(node))
# function join_multiblock(node::Node)
#     @assert node.element isa Documenter.MultiCodeBlock
#     io = IOBuffer()
#     codeblocks = [n.element::MarkdownAST.CodeBlock for n in node.children]
#     for (i, thing) in enumerate(codeblocks)
#         print(io, thing.code)
#         if i != length(codeblocks)
#             println(io)
#             if findnext(x -> x.info == node.element.language, codeblocks, i + 1) == i + 1
#                 println(io)
#             end
#         end
#     end
#     return MarkdownAST.CodeBlock(node.element.language, String(take!(io)))
# end

# function _print_code_escapes_minted(io, s::AbstractString)
#     for ch in s
#         ch === '#' ? _print(io, "##%") :
#         ch === '%' ? _print(io, "#%%") : # Note: "#\\%%" results in pygmentize error...
#         ch === '⊻' ? _print(io, "#\\unicodeveebar%") :
#         _print(io, ch)
#     end
#     return
# end

function typst(io::Context, node::Node, ::MarkdownAST.Paragraph)
    typst(io, node.children)
    _println(io, "\n")
    return nothing
end

# function typst(io::Context, node::Node, ::MarkdownAST.BlockQuote)
#     wrapblock(io, "quote") do
#         typst(io, node.children)
#     end
#     return
# end

function typst(io::Context, node::Node, md::MarkdownAST.Admonition)
    if md.category ∉ ("note", "tip", "important", "warning", "caution")
        warn("unsupported admonition type $(md.category)")
    else
        _println(io, "#$(md.category)[")
        _print(io, "    ")
        typst(io, node.children)
        _println(io, "]")
    end
    return nothing
end

# function typst(io::Context, node::Node, f::MarkdownAST.FootnoteDefinition)
#     id = get(io.footnotes, f.id, 1)
#     _print(io, "\\footnotetext[", id, "]{")
#     typst(io, node.children)
#     _println(io, "}")
#     return
# end

# function typst(io::Context, node::Node, list::MarkdownAST.List)
#     # TODO: MarkdownAST doesn't support lists starting at arbitrary numbers
#     isordered = (list.type === :ordered)
#     ordered = (list.type === :bullet) ? -1 : 1
#     # `\begin{itemize}` is used here for both ordered and unordered lists since providing
#     # custom starting numbers for enumerated lists is simpler to do by manually assigning
#     # each number to `\item` ourselves rather than using `\setcounter{enumi}{<start>}`.
#     #
#     # For an ordered list starting at 5 the following will be generated:
#     #
#     # \begin{itemize}
#     #   \item[5. ] ...
#     #   \item[6. ] ...
#     #   ...
#     # \end{itemize}
#     #
#     pad = ndigits(ordered + length(node.children)) + 2
#     fmt = n -> (isordered ? "[$(rpad("$(n + ordered - 1).", pad))]" : "")
#     wrapblock(io, "itemize") do
#         for (n, item) in enumerate(node.children)
#             _print(io, "\\item$(fmt(n)) ")
#             typst(io, item.children)
#             n < length(node.children) && _println(io)
#         end
#     end
#     return
# end

# function typst(io::Context, node::Node, e::MarkdownAST.ThematicBreak)
#     _println(io, "{\\rule{\\textwidth}{1pt}}")
#     return
# end

# This (equation*, split) math env seems to be the only way to correctly render all the
# equations in the Julia manual. However, if the equation is already wrapped in
# align/align*, then there is no need to further wrap it (in fact, it will break).
# function typst(io::Context, node::Node, math::MarkdownAST.DisplayMath)
#     if occursin(r"^\\begin\{align\*?\}", math.math)
#         _print(io, math.math)
#     else
#         _print(io, "\\begin{equation*}\n\\begin{split}")
#         _print(io, math.math)
#         _println(io, "\\end{split}\\end{equation*}")
#     end
#     return
# end

# function typst(io::Context, node::Node, table::MarkdownAST.Table)
#     rows = MarkdownAST.tablerows(node)
#     _println(io, "\n\\begin{table}[h]\n\\centering")
#     _print(io, "\\begin{tabulary}{\\linewidth}")
#     _println(io, "{", uppercase(join(spec_to_align.(table.spec), ' ')), "}")
#     _println(io, "\\toprule")
#     for (i, row) in enumerate(rows)
#         for (j, cell) in enumerate(row.children)
#             j === 1 || _print(io, " & ")
#             typst(io, cell.children)
#         end
#         _println(io, " \\\\")
#         if i === 1
#             _println(io, "\\toprule")
#         end
#     end
#     _println(io, "\\bottomrule")
#     _println(io, "\\end{tabulary}\n")
#     _println(io, "\\end{table}\n")
#     return
# end
# spec_to_align(spec::Symbol) = Symbol(first(String(spec)))

# function typst(io::Context, node::Node, raw::Documenter.RawNode)
#     raw.name === :latex && _println(io, "\n", raw.text, "\n")
#     return
# end

function typst(io::Context, node::Node, e::MarkdownAST.Text)
    return latexesc(io, e.text)
end

# function typst(io::Context, node::Node, image::Documenter.LocalImage)
#     # TODO: also print the .title field somehow
#     wrapblock(io, "figure") do
#         _println(io, "\\centering")
#         wrapinline(io, "includegraphics[max width=\\linewidth]") do
#             _print(io, replace(image.path, "\\" => "/"))
#         end
#         _println(io)
#         wrapinline(io, "caption") do
#             typst(io, node.children)
#         end
#         _println(io)
#     end
#     return
# end

# function typst(io::Context, node::Node, image::MarkdownAST.Image)
#     # TODO: also print the .title field somehow
#     wrapblock(io, "figure") do
#         _println(io, "\\centering")
#         @warn "images with absolute URLs not supported in LaTeX output in $(Documenter.locrepr(io.filename))" url = image.destination
#         # We nevertheless output an \includegraphics with the URL. The LaTeX build will
#         # then give an error, indicating to the user that something wrong.
#         url = replace(image.destination, "\\" => "/") # use / on Windows too.
#         wrapinline(io, "includegraphics[max width=\\linewidth]") do
#             _print(io, url)
#         end
#         _println(io)
#         wrapinline(io, "caption") do
#             typst(io, node.children)
#         end
#         _println(io)
#     end
#     return
# end

# function typst(io::Context, node::Node, f::MarkdownAST.FootnoteLink)
#     id = get!(io.footnotes, f.id, length(io.footnotes) + 1)
#     _print(io, "\\footnotemark[", id, "]")
#     return
# end

# function typst(io::Context, node::Node, link::Documenter.PageLink)
#     # If we're in a header, we don't want to print any \hyperlinkref commands,
#     # so we handle this here.
#     if io.in_header
#         typst(io, node.children)
#         return
#     end
#     # This branch is the normal case, when we're not in a header.
#     # TODO: this link handling does not seem correct
#     if !isempty(link.fragment)
#         id = _hash(link.fragment)
#         wrapinline(io, "hyperlinkref") do
#             _print(io, id)
#         end
#     else
#         wrapinline(io, "href") do
#             path = Documenter.pagekey(io.doc, link.page)
#             latexesc(io, path)
#         end
#     end
#     _print(io, "{")
#     typst(io, node.children)
#     _print(io, "}")
#     return
# end

# function typst(io::Context, node::Node, link::Documenter.LocalLink)
#     # If we're in a header, we don't want to print any \hyperlinkref commands,
#     # so we handle this here.
#     if io.in_header
#         typst(io, node.children)
#         return
#     end
#     # This branch is the normal case, when we're not in a header.
#     # TODO: this link handling does not seem correct
#     wrapinline(io, "href") do
#         href = isempty(link.fragment) ? link.path : "$(link.path)#($(link.fragment))"
#         latexesc(io, href)
#     end
#     _print(io, "{")
#     typst(io, node.children)
#     _print(io, "}")
#     return
# end

# function typst(io::Context, node::Node, link::MarkdownAST.Link)
#     # If we're in a header, we don't want to print any \hyperlinkref commands,
#     # so we handle this here.
#     if io.in_header
#         typst(io, node.children)
#         return
#     end
#     # This branch is the normal case, when we're not in a header.
#     # TODO: handle the .title attribute
#     wrapinline(io, "href") do
#         latexesc(io, link.destination)
#     end
#     _print(io, "{")
#     typst(io, node.children)
#     _print(io, "}")
#     return
# end

# Metadata Nodes get dropped from the final output for every format but are needed throughout
# rest of the build and so we just leave them in place and print a blank line in their place.
# typst(io::Context, node::Node, ::Documenter.MetaNode) = _println(io, "\n")

# # In the original AST, SetupNodes were just mapped to empty Markdown.MD() objects.
# typst(io::Context, node::Node, ::Documenter.SetupNode) = nothing

# function typst(io::Context, node::Node, value::MarkdownAST.JuliaValue)
#     @warn(
#         """
#         Unexpected Julia interpolation in the Markdown. This probably means that you
#         have an unbalanced or un-escaped \$ in the text.

#         To write the dollar sign, escape it with `\\\$`

#         We don't have the file or line number available, but we got given the value:

#         `$(value.ref)` which is of type `$(typeof(value.ref))`
#         """
#     )
#     return latexesc(io, string(value.ref))
# end

# TODO: Implement SoftBreak, Backslash (but they don't appear in standard library Markdown conversions)
# typst(io::Context, node::Node, ::MarkdownAST.LineBreak) = _println(io, "\\\\")

# Documenter.

const typst_escape_characters = Dict{Char,AbstractString}(
    '#' => "\\#", '@' => "\\@", '$' => "\\\$", '[' => "\\[", ']' => "\\]"
)

latexesc(io, ch::AbstractChar) = _print(io, get(typst_escape_characters, ch, ch))

function latexesc(io, s::AbstractString)
    for ch in s
        latexesc(io, ch)
    end
    return nothing
end

latexesc(s) = sprint(latexesc, s)

function files!(out::Vector, v::Vector, depth)
    for each in v
        files!(out, each, depth + 1)
    end
    return out
end

# Tuples come from `hide(page)` with either
# (visible, nothing,    page,         children) or
# (visible, page.first, pages.second, children)
function files!(out::Vector, v::Tuple, depth)
    files!(out, v[2] == nothing ? v[3] : v[2] => v[3], depth)
    files!(out, v[4], depth)
    return nothing
end

files!(out, s::AbstractString, depth) = push!(out, ("", s, depth))

function files!(out, p::Pair{<:AbstractString,<:Any}, depth)
    # Hack time. Because of Julia's typing, something like
    # `"Introduction" => "index.md"` may get typed as a `Pair{String,Any}`!
    if p[2] isa AbstractString
        push!(out, (p.first, p.second, depth))
    else
        push!(out, (p.first, "", depth))
        files!(out, p.second, depth)
    end
    return out
end

files(v::Vector) = files!(Tuple{String,String,Int}[], v, 0)

end
