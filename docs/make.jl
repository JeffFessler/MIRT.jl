using Documenter, DemoCards
using MIRT

# examples_templates, examples_theme = cardtheme("grid")
# examples, examples_cb = makedemos("examples", examples_templates)

format = Documenter.HTML(
    edit_link = "master",
    prettyurls = get(ENV, "CI", nothing) == "true",
#   assets = [examples_theme],
)


makedocs(
    modules  = [MIRT],
    format   = format,
    sitename = "MIRT.jl",
    pages    = [
        "Home" => "index.md",
        "Table of Contents" => "toc.md",
        # "Examples" => examples,
        "Function References" => "reference.md",
    ],
    authors = "Jeff Fessler and contributors",
)

# examples_cb()

deploydocs(
    repo = "github.com/JeffFessler/MIRT.jl.git",
)
