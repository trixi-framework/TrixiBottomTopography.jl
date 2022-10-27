# Include packages
using TrixiBottomTopography
using Test

# Start testing
@testset "TrixiBottomTopography.jl" begin
    
  # Test data conversion
  include("test_convert.jl")

  # Test one dimensional B-spline interpolation
  # include("test_1d_b_spline.jl")

  # Test two dimensional B-spline interpolation
  # include("test_2d_b_spline.jl")

end
