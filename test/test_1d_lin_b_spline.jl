using TrixiBottomTopography

# Define data path
data = "data/rhine_1d_10_x_1.txt"

# Define B-spline structure
spline_struct = linear_b_spline(data)
# Define B-spline interpolation function
spline_func(x) = spline_interpolation(spline_struct, x)

# Define interpolation points
n = 100
x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))

# Get interpolated values
y_int_pts = spline_func.(x_int_pts)