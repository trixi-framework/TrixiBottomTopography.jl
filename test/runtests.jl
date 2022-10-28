# Include packages
using TrixiBottomTopography
using Test

# Start testing
@testset "TrixiBottomTopography.jl" begin
    
  # Test data conversion
  @testset "Conversion" begin
    include("test_convert.jl")
  end
  
  @testset "Linear B-spline interpolation" begin
    include("test_1d_lin_b_spline.jl")
  end
  
  # Test two dimensional B-spline interpolation
  # include("test_2d_b_spline.jl")

end
