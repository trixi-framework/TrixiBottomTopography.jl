using Documenter
using TrixiBottomTopography

makedocs(
    sitename = "TrixiBottomTopography.jl",
    format = Documenter.HTML(),
    modules = [TrixiBottomTopography],
    pages = [
        "Home" => "index.md",
        "API" => "reference.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
