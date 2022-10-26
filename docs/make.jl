using Documenter
import Pkg
using TrixiBottomTopography

# Define module-wide setups such that the respective modules are available in doctests
DocMeta.setdocmeta!(TrixiBottomTopography,     :DocTestSetup, :(using TrixiBottomTopography);     recursive=true)

# Make documentation
makedocs(
    # Specify modules for which docstrings should be shown
    modules = [TrixiBottomTopography],
    # Set sitename to TrixiBottomTopography
    sitename="TrixiBottomTopography.jl",
    # Provide additional formatting options
    format = Documenter.HTML(
        # Disable pretty URLs during manual testing
        prettyurls = get(ENV, "CI", nothing) == "true",
        # Explicitly add favicon as asset
        # assets = ["assets/favicon.ico"],
        # Set canonical URL to GitHub pages URL
        canonical = "https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev"
    ),
    # Explicitly specify documentation structure
    pages = [
        "Home" => "index.md",
        "Reference" => "reference.md",
        "License" => "license.md"
    ],
    # strict = true # to make the GitHub action fail when doctests fail
    strict = Documenter.except(:cross_references)
)

deploydocs(
    repo = "github.com/maxbertrand1996/TrixiBottomTopography.jl",
    devbranch = "main",
    push_preview = true
)