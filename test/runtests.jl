# Include packages
using TrixiBottomTopography
using Test

# Start testing
@testset "TrixiBottomTopography.jl" begin
    
  # Test data conversion
  @testset "Conversion" begin
    include("test_convert.jl")
  end
  
  # Linear B-splines
  @testset "Linear B-spline interpolation" begin
    include("test_linear_b_spline.jl")
  end
  
  # Cubic B-splines
  @testset "Cubic B-spline interpolation" begin
    include("test_cubic_b_spline.jl")
  end

  # Bilinear B-splines
  @testset "Bilinear B-spline interpolation" begin
    include("test_bilinear_b_spline.jl")
  end

  # Bicubic B-splines
  @testset "Bicubic B-spline interpolation" begin
    include("test_bicubic_b_spline.jl")
  end

end
