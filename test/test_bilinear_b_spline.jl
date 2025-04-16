using TrixiBottomTopography

# Define data path
data = joinpath(@__DIR__, "data", "rhine_2d_10.txt")

# Define B-spline structure
spline_struct = BilinearBSpline(data)
# Define B-spline interpolation function
spline_func(x, y) = spline_interpolation(spline_struct, x, y)

# Test function at arbitrary point
@test 39.26 <= spline_func(357555, 5646555) <= 39.26000000000001
