##############################################################################
# Script which uses the functionalities implemented in TrixiBottomTopography #
# to plot a one dimensional cubic B-spline interpolated section of the Rhine #
# river with not-a-knot end condition, smoothing with data extracted along   #
# the x-direction.                                                           #
##############################################################################

# Include packages
using TrixiBottomTopography

# Define data path
root_dir = pkgdir(TrixiBottomTopography)
data = joinpath(root_dir, "examples", "data", "rhine_data_1d_20_x.txt")

# Define B-spline structure
spline_struct = CubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 999)
# Define B-spline interpolation function
spline_func(x) = spline_interpolation(spline_struct, x)

# Evaluate the cubic B-spline on a new set of nodes and plot
if isdefined(Main, :Makie)
  # Define interpolation points
  n = 200
  x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))

  # Get interpolated values
  y_int_pts = spline_func.(x_int_pts)

  # Get the original interpolation knots
  x_knots = spline_struct.x
  y_knots = spline_func.(x_knots)

  plot_topography_with_interpolation_knots(x_int_pts, y_int_pts, x_knots, y_knots;
                                           xlabel = "ETRS89 East", ylabel = "DHHN2016 Height" )
end
