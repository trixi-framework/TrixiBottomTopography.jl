using TrixiBottomTopography

# Define data path
data = "data/rhine_1d_10_x_1.txt"

# Define B-spline structure
spline_struct = linear_b_spline(data)
# Define B-spline interpolation function
spline_func(x) = spline_interpolation(spline_struct, x)

# Test function at arbitrary point
@test spline_func(357555) == 46.19