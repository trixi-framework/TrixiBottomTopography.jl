using TrixiBottomTopography

# Free end condition
@testset "Free end" begin

  # Define data path
  root_dir = pkgdir(TrixiBottomTopography)
  data = joinpath(root_dir, "test", "data", "rhine_1d_10_x_1.txt")

  # Define B-spline structure
  spline_struct = CubicBSpline(data)
  # Define B-spline interpolation function
  spline_func(x) = spline_interpolation(spline_struct, x)

  # Test function at arbitrary point
  @test 46.19366216271 < spline_func(357555) < 46.19366216272

end


# Free end condition with smoothing
@testset "Free end + smoothing" begin

  # Define data path
  root_dir = pkgdir(TrixiBottomTopography)
  data = joinpath(root_dir, "test", "data", "rhine_1d_10_x_100.txt")

  # Define B-spline structure
  spline_struct = CubicBSpline(data; smoothing_factor = 9999)
  # Define B-spline interpolation function
  spline_func(x) = spline_interpolation(spline_struct, x)

  # Test function at arbitrary point
  @test 47.85931152 < spline_func(357555) < 47.85931153

end

# Not-a-knot end condition
@testset "Not-a-knot end" begin

  # Define data path
  root_dir = pkgdir(TrixiBottomTopography)
  data = joinpath(root_dir, "test", "data", "rhine_1d_10_y_1.txt")

  # Define B-spline structure
  spline_struct = CubicBSpline(data; end_condition = "not-a-knot")
  # Define B-spline interpolation function
  spline_func(x) = spline_interpolation(spline_struct, x)

  # Test function at arbitrary point
  @test 47.92937967318 < spline_func(5646555) < 47.92937967319

end

# Not-a-knot end condition with smoothing
@testset "Not-a-knot end + smoothing" begin

  # Define data path
  root_dir = pkgdir(TrixiBottomTopography)
  data = joinpath(root_dir, "test", "data", "rhine_1d_10_y_100.txt")

  # Define B-spline structure
  spline_struct = CubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 9999)
  # Define B-spline interpolation function
  spline_func(x) = spline_interpolation(spline_struct, x)

  # Test function at arbitrary point
  @test 47.71995039 < spline_func(5646555) < 47.71995040

end
