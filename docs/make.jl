using Documenter
using RMPBT
DocMeta.setdocmeta!(RMPBT, :DocTestSetup, :(using RMPBT); recursive=true)
makedocs(
    sitename = "RMPBT.jl",
    format = Documenter.HTML(),
    modules = [RMPBT],
    authors="Sai Suraj Karra <saisurajkarra.com> and contributors",
    repo="https://github.com//RMPBT.jl/blob/{commit}{path}#{line}",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://saisurajkarra.github.io/RMPBT.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
deploydocs(;
    repo="github.com/saisurajkarra/RMPBT.jl",
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#

