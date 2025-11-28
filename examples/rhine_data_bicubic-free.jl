##############################################################################
# Script which uses the functionalities implemented in TrixiBottomTopography #
# to plot a bicubic B-spline interpolated section of the Rhine river with    #
# free end condition.                                                        #
##############################################################################

# Include packages
using TrixiBottomTopography

# Define data path
root_dir = pkgdir(TrixiBottomTopography)
data = joinpath(root_dir, "examples", "data", "rhine_data_2d_20.txt")

# Define B-spline structure
spline_struct = BicubicBSpline(data)
# Define B-spline interpolation function
spline_func(x, y) = spline_interpolation(spline_struct, x, y)

# Evaluate the bicubic B-spline on a new set of nodes and plot
if isdefined(Main, :Makie)
    # Define interpolation points
    n = 100
    x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))
    y_int_pts = Vector(LinRange(spline_struct.y[1], spline_struct.y[end], n))

    # Get interpolated matrix
    z_int_pts = evaluate_bicubicspline_interpolant(spline_func, x_int_pts, y_int_pts)

    # Get the original interpolation knots
    x_knots = spline_struct.x
    y_knots = spline_struct.y
    z_knots = evaluate_bicubicspline_interpolant(spline_func, x_knots, y_knots)

    plot_topography_with_interpolation_knots(x_int_pts, y_int_pts, z_int_pts,
                                             x_knots, y_knots, z_knots;
                                             xlabel = "ETRS89 East",
                                             ylabel = "DHHN2016 Height", zlabel = "H")
end