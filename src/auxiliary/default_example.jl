"""
    default_example()

Function which calls the example file "rhine_data_bicubic-nak.jl" as a default example to
check if including the package has worked creates a bicubic B-spline
with the "not-a-knot" boundary condition.
If a [Makie.jl](https://github.com/JuliaPlots/Makie.jl/) backend
such as [GLMakie.jl](https://github.com/JuliaPlots/GLMakie.jl/)
is included a plot will be generated as well.
"""
function default_example()

  # Location of example file
  root_dir = pkgdir(TrixiBottomTopography)
  example_file = joinpath(root_dir, "examples", "rhine_data_bicubic-nak.jl")

  # Include example
  include(example_file)

end
