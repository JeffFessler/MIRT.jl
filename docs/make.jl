
using MIRT
using Documenter
#using DemoCards
using Literate

# examples_templates, examples_theme = cardtheme("grid")
# examples, examples_cb = makedemos("examples", examples_templates)

# generate examples using Literate
lit = joinpath(@__DIR__, "lit")
src = joinpath(@__DIR__, "src")
notebooks = joinpath(src, "notebooks")

ENV["GKS_ENCODING"] = "utf-8"

DocMeta.setdocmeta!(MIRT, :DocTestSetup, :(using MIRT); recursive=true)

execute = true # Set to true for executing notebooks and documenter
nb = false # Set to true to generate the notebooks
for (root, _, files) in walkdir(lit), file in files
    splitext(file)[2] == ".jl" || continue
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, lit=>src))[1]
    Literate.markdown(ipath, opath, documenter = execute)
    nb && Literate.notebook(ipath, notebooks, execute = execute)
end

# Documentation structure
ismd(f) = splitext(f)[2] == ".md"
pages(folder) =
    [joinpath(folder, f) for f in readdir(joinpath(src, folder)) if ismd(f)]

@show isci = get(ENV, "CI", nothing) == "true"

format = Documenter.HTML(;
    prettyurls = isci,
#   edit_link = "master",
#   assets = [examples_theme],
)

makedocs(;
    modules = [MIRT],
    authors = "Jeff Fessler and contributors",
    sitename = "MIRT.jl",
    format,
    pages = [
        "Home" => "index.md",
        "Table of Contents" => "toc.md",
        # "Examples" => examples,
        "Examples" => pages("examples"),
        "Function References" => "reference.md",
    ],
)

# examples_cb()

if isci
    deploydocs(;
        repo = "github.com/JeffFessler/MIRT.jl.git",
        devbranch = "master",
        versions = ["stable" => "v^", "dev" => "dev"],
    )
else
    @warn "may need to rm -r src/examples"
end
