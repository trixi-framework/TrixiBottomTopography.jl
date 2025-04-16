##############################################################################
# Script which uses the functionalities implemented in TrixiBottomTopography #
# to plot a one dimensional cubic B-spline interpolated section of the Rhine #
# river with not-a-knot end condition.                                       #
##############################################################################

# Include packages
using TrixiBottomTopography
using Plots

# Define data path
data = joinpath(@__DIR__, "data", "rhine_data_1d_20_x.txt")

# Define B-spline structure
spline_struct = CubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 999)
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
     title="Cubic B-spline interpolation\nwith not-a-knot end condition and smoothing")