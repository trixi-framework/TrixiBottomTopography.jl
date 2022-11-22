using Documenter
import Pkg
using TrixiBottomTopography

# Define module-wide setups such that the respective modules are available in doctests
DocMeta.setdocmeta!(TrixiBottomTopography, :DocTestSetup, :(using TrixiBottomTopography);
                    recursive=true)

# Make documentation
makedocs(
    sitename = "TrixiBottomTopography.jl",
    format = Documenter.HTML(
        # Set canonical URL to GitHub pages URL
        canonical = "https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev"
    ),
    modules = [TrixiBottomTopography],
    pages = [
        "Home" => "index.md",
        "Overview" => [
            "Data Conversion" => "conversion.md",
            "B-spline structure" => "structure.md",
            "B-spline function" => "function.md",
        ],
        "Trixi.jl examples" => [
            "1D dam break" => "Trixi_1D.md",
            # "2D dam break" => "Trixi_2D.md",
        ],
        "Reference" => "reference.md",
        "Licence" => "licence.md"
    ],
    # strict = true # to make the GitHub action fail when doctests fail
    strict = Documenter.except(:cross_references)
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://github.com/maxbertrand1996/TrixiBottomTopography.jl",
    devbranch = "main",
    push_preview = true
)
