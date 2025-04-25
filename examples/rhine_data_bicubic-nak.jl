##############################################################################
# Script which uses the functionalities implemented in TrixiBottomTopography #
# to plot a bicubic B-spline interpolated section of the Rhine river with    #
# not-a-knot end condition and smoothing.                                    #
##############################################################################

# Include packages
using TrixiBottomTopography

# Define data path
root_dir = pkgdir(TrixiBottomTopography)
data = joinpath(root_dir, "examples", "data", "rhine_data_2d_20.txt")

# Define B-spline structure
spline_struct = BicubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 9999)
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

    plot_topography(x_int_pts, y_int_pts, z_int_pts;
                    xlabel = "ETRS89\n East",
                    ylabel = "ETRS89\n North",
                    zlabel = "DHHN2016\n Height",
                    azimuth_angle = 54 * pi / 180,
                    elevation_angle = 27 * pi / 180)
end
