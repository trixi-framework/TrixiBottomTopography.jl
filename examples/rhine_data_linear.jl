###############################################################################
# Script which uses the functionalities implemented in TrixiBottomTopography  #
# to polt a one dimensional linear B-spline interpolated section of the Rhine #
# river with not-a-knot end condition.                                        #
###############################################################################

# Include packages
using TrixiBottomTopography
using Plots

# Get root directory
dir_path = pkgdir(TrixiBottomTopography)

# Define data path
data = string(dir_path, "/examples/data/rhine_data_1d_20_y.txt")

# Define B-spline structure
spline_struct = LinearBSpline(data)
# Define B-spline interpolation function
spline_func(x) = spline_interpolation(spline_struct, x)

# Define interpolation points
n = 100
x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))

# Get interpolated values
y_int_pts = spline_func.(x_int_pts)

# Plotting
pyplot()
plot(x_int_pts, y_int_pts,
     xlabel="ETRS89 East", ylabel="DHHN2016 Height",
     label="Bottom topography", 
     title="Linear B-spline interpolation")