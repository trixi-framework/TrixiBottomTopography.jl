using Documenter
import Pkg
using TrixiBottomTopography

# Copy list of authors to not need to synchronize it manually.
# Since the authors header exists twice we create a unique identifier for the docs section.
authors_text = read(joinpath(dirname(@__DIR__), "AUTHORS.md"), String)
authors_text = replace(authors_text,
                       "in the [LICENSE.md](LICENSE.md) file" => "under [License](@ref)",
                       "# Authors" => "# [Authors](@id trixi_bt_authors)")
write(joinpath(@__DIR__, "src", "authors.md"), authors_text)

# Define module-wide setups such that the respective modules are available in doctests
DocMeta.setdocmeta!(TrixiBottomTopography, :DocTestSetup, :(using TrixiBottomTopography);
                    recursive=true)

# Copy code of conduct to not need to synchronize it manually
code_of_conduct_text = read(joinpath(dirname(@__DIR__), "CODE_OF_CONDUCT.md"), String)
code_of_conduct_text = replace(code_of_conduct_text,
                               "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref trixi_bt_authors)")
write(joinpath(@__DIR__, "src", "code_of_conduct.md"), code_of_conduct_text)


# Copy contributing information to not need to synchronize it manually
contributing_text = read(joinpath(dirname(@__DIR__), "CONTRIBUTING.md"), String)
contributing_text = replace(contributing_text,
                            "[LICENSE.md](LICENSE.md)" => "[License](@ref)",
                            "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref trixi_bt_authors)")
write(joinpath(@__DIR__, "src", "contributing.md"), contributing_text)

# Copy contents form README to the starting page to not need to synchronize it manually
readme_text = read(joinpath(dirname(@__DIR__), "README.md"), String)
readme_text = replace(readme_text,
                      "[LICENSE.md](LICENSE.md)" => "[License](@ref)",
                      "[AUTHORS.md](AUTHORS.md)" => "[Authors](@ref trixi_bt_authors)",
                      "<p" => "```@raw html\n<p",
                      "p>" => "p>\n```",
                      r"\[comment\].*\n" => "")    # remove comments
write(joinpath(@__DIR__, "src", "home.md"), readme_text)

# Make documentation
makedocs(;
         modules = [TrixiBottomTopography],
         authors = "Andrew R. Winters <andrew.ross.winters@liu.se>, Michael Schlottke-Lakemper <michael@sloede.com>",
         repo = "https://github.com/trixi-framework/TrixiBottomTopography.jl/blob/{commit}{path}#{line}",
         sitename = "TrixiBottomTopography.jl",
         format = Documenter.HTML(;
                                  # Disable pretty URLs during manual testing
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  # Set canonical URL to GitHub pages URL
                                  canonical = "https://trixi-framework.github.io/TrixiBottomTopography.jl/stable",
                                  edit_link = "main",
                                  size_threshold_ignore = ["index.md"],),
         # Explicitly specify documentation structure
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
         "Licence" => "licence.md"],
        )

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/trixi-framework/TrixiBottomTopography.jl.git",
    devbranch = "main",
    push_preview = true
)
