using Documenter
using RMPBT

makedocs(
    sitename = "RMPBT",
    format = Documenter.HTML(),
    modules = [RMPBT]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
