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
        canonical = "https://trixi-framework.github.io/TrixiBottomTopography.jl/dev"
    ),
    modules = [TrixiBottomTopography],
    pages = [
        "Home" => "index.md",
        "Overview" => [
            "Data conversion" => "conversion.md",
            "B-spline structure" => "structure.md",
            "B-spline function" => "function.md",
        ],
        "Trixi.jl examples" => "Trixi.md",
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
    repo = "https://trixi-framework.github.io/TrixiBottomTopography.jl",
    devbranch = "main",
    push_preview = true
)
