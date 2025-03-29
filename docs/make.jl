using Documenter
import Pkg
using TrixiBottomTopography

# Copy list of authors to not need to synchronize it manually
authors_text = read(joinpath(dirname(@__DIR__), "AUTHORS.md"), String)
authors_text = replace(authors_text,
                       "in the [LICENSE.md](LICENSE.md) file" => "under [License](@ref)")
write(joinpath(@__DIR__, "src", "authors.md"), authors_text)

# Define module-wide setups such that the respective modules are available in doctests
DocMeta.setdocmeta!(TrixiBottomTopography, :DocTestSetup, :(using TrixiBottomTopography);
                    recursive=true)

# Copy some files from the repository root directory to the docs and modify them
# as necessary
# Based on: https://github.com/ranocha/SummationByPartsOperators.jl/blob/0206a74140d5c6eb9921ca5021cb7bf2da1a306d/docs/make.jl#L27-L41
open(joinpath(@__DIR__, "src", "code_of_conduct.md"), "w") do io
    # Point to source license file
    println(io,
            """
            ```@meta
            EditURL = "https://github.com/trixi-framework/TrixiBottomTopography.jl/blob/main/CODE_OF_CONDUCT.md"
            ```
            """)
    # Write the modified contents
    println(io, "# [Code of Conduct](@id code-of-conduct)")
    println(io, "")
    for line in eachline(joinpath(dirname(@__DIR__), "CODE_OF_CONDUCT.md"))
        line = replace(line, "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)")
        println(io, "> ", line)
    end
end

open(joinpath(@__DIR__, "src", "contributing.md"), "w") do io
    # Point to source license file
    println(io,
            """
            ```@meta
            EditURL = "https://github.com/trixi-framework/TrixiBottomTopography.jl/blob/main/CONTRIBUTING.md"
            ```
            """)
    # Write the modified contents
    for line in eachline(joinpath(dirname(@__DIR__), "CONTRIBUTING.md"))
        line = replace(line, "[LICENSE.md](LICENSE.md)" => "[License](@ref)")
        line = replace(line, "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref)")
        println(io, line)
    end
end

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
        "Authors" => "authors.md",
        "Contributing" => "contributing.md",
        "Code of Conduct" => "code_of_conduct.md",
        "Licence" => "licence.md"
    ],
    strict = true # to make the GitHub action fail when doctests fail
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/trixi-framework/TrixiBottomTopography.jl.git",
    devbranch = "main",
    push_preview = true
)
