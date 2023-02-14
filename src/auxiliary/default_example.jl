"""
    TBT_default_example()

Function which calls the example file "rhine_data_bicubic-nak.jl" as a default example to
check if including the package has worked and to see a quick example.
"""
function TBT_default_example()
  
  # Get root directory
  dir_path = pkgdir(TrixiBottomTopography)

  # Location of example file
  example_file = string(dir_path, "/examples/rhine_data_bicubic-nak.jl")

  # Include example
  include(example_file)

end
