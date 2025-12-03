using Documenter
import Pkg
using TrixiBottomTopography
using Changelog: Changelog
# To load the extensions
import CairoMakie
import KernelInterpolation

# Copy list of authors to not need to synchronize it manually.
# Since the authors header exists twice we create a unique identifier for the docs section.
authors_text = read(joinpath(dirname(@__DIR__), "AUTHORS.md"), String)
authors_text = replace(authors_text,
                       "in the [LICENSE.md](LICENSE.md) file" => "under [License](@ref)")
authors_text = replace(authors_text,
                       "# Authors" => "# [Authors](@id trixi_bt_authors)")
write(joinpath(@__DIR__, "src", "authors.md"), authors_text)

# Define module-wide setups such that the respective modules are available in doctests
DocMeta.setdocmeta!(TrixiBottomTopography, :DocTestSetup, :(using TrixiBottomTopography);
                    recursive = true)

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
        line = replace(line,
                       "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref trixi_bt_authors)")
        println(io, "> ", line)
    end
end

# Copy contributing information to not need to synchronize it manually
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
        line = replace(line,
                       "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref trixi_bt_authors)")
        println(io, line)
    end
end

# Copy contents form README to the starting page to not need to synchronize it manually
readme_text = read(joinpath(dirname(@__DIR__), "README.md"), String)
readme_text = replace(readme_text,
                      "[LICENSE.md](LICENSE.md)" => "[License](@ref)")
readme_text = replace(readme_text,
                      "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref trixi_bt_authors)")
readme_text = replace(readme_text,
                      "<p" => "```@raw html\n<p")
readme_text = replace(readme_text,
                      "p>" => "p>\n```")
readme_text = replace(readme_text,
                      r"\[comment\].*\n" => "")    # remove comments
write(joinpath(@__DIR__, "src", "index.md"), readme_text)

# Create changelog
Changelog.generate(Changelog.Documenter(),                            # output type
                   joinpath(@__DIR__, "..", "NEWS.md"),               # input file
                   joinpath(@__DIR__, "src", "changelog_tmp.md");     # output file
                   repo = "trixi-framework/TrixiBottomTopography.jl", # default repository for links
                   branch = "main",)
# Fix edit URL of changelog
open(joinpath(@__DIR__, "src", "changelog.md"), "w") do io
    for line in eachline(joinpath(@__DIR__, "src", "changelog_tmp.md"))
        if startswith(line, "EditURL")
            line = "EditURL = \"https://github.com/trixi-framework/TrixiBottomTopography.jl/blob/main/NEWS.md\""
        end
        println(io, line)
    end
end
# Remove temporary file
rm(joinpath(@__DIR__, "src", "changelog_tmp.md"))

# Make documentation
makedocs(;
         modules = [TrixiBottomTopography,
             Base.get_extension(TrixiBottomTopography, :MakieExt),
             Base.get_extension(TrixiBottomTopography, :KernelInterpolationExt)],
         authors = "Andrew R. Winters <andrew.ross.winters@liu.se>, Michael Schlottke-Lakemper <michael@sloede.com>",
         sitename = "TrixiBottomTopography.jl",
         format = Documenter.HTML(;
                                  # Disable pretty URLs during manual testing
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  # Set canonical URL to GitHub pages URL
                                  canonical = "https://trixi-framework.github.io/TrixiBottomTopography.jl/stable",
                                  edit_link = "main",
                                  size_threshold_ignore = ["index.md"],),
         # Explicitly specify documentation structure
         pages = ["Home" => "index.md",
             "Overview" => [
                 "Data conversion" => "conversion.md",
                 "B-spline structure" => "structure.md",
                 "B-spline function" => "function.md"
             ],
             "TrixiShallowWater.jl examples" => "trixishallowwater_jl_examples.md",
             "Advanced topics & developers" => ["Development" => "development.md",
                 "Style guide" => "styleguide.md",
                 "Testing" => "testing.md"],
             "Reference" => "reference.md",
             "Changelog" => "changelog.md",
             "Authors" => "authors.md",
             "Contributing" => "contributing.md",
             "Code of Conduct" => "code_of_conduct.md",
             "License" => "license.md"])

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(repo = "github.com/trixi-framework/TrixiBottomTopography.jl.git",
           devbranch = "main",
           # Only push previews if all the relevant environment variables are non-empty.
           push_preview = all(!isempty,
                              (get(ENV, "GITHUB_TOKEN", ""),
                               get(ENV, "DOCUMENTER_KEY", ""))))
