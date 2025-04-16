"""
    TBT_default_example()

Function which calls the example file "rhine_data_bicubic-nak.jl" as a default example to
check if including the package has worked and to see a quick example.
"""
function TBT_default_example()

  # Location of example file
  example_file = joinpath(@__DIR__, "examples", "rhine_data_bicubic-nak.jl")

  # Include example
  include(example_file)

end
