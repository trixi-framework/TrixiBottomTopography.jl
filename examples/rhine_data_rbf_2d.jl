##############################################################################
# Script which uses the functionalities implemented in TrixiBottomTopography #
# to plot a RBF interpolated section of the Rhine river.                     #
##############################################################################

# Include packages
using TrixiBottomTopography
using KernelInterpolation

# Define data path
root_dir = pkgdir(TrixiBottomTopography)
data = joinpath(root_dir, "examples", "data", "rhine_data_2d_20.txt")

# Define RBF structure
kernel = ThinPlateSplineKernel{2}()
itp = RBFInterpolation2D(data, kernel)
# Define interpolation function
itp_func(x, y) = itp([x; y])

# Evaluate the thin plate spline on a new set of nodes and plot
if isdefined(Main, :Makie)
    nodes = nodeset(itp)
    x_knots = sort(unique(nodes[:, 1]))
    y_knots = sort(unique(nodes[:, 2]))
    z_knots = evaluate_two_dimensional_interpolant(itp_func, x_knots, y_knots)

    # Define interpolation points
    n = 100
    x_int_pts = Vector(LinRange(x_knots[1], x_knots[end], n))
    y_int_pts = Vector(LinRange(y_knots[1], y_knots[end], n))

    # Get interpolated matrix
    z_int_pts = evaluate_two_dimensional_interpolant(itp_func, x_int_pts, y_int_pts)

    plot_topography_with_interpolation_knots(x_int_pts, y_int_pts, z_int_pts,
                                             x_knots, y_knots, z_knots;
                                             xlabel = "ETRS89 East",
                                             ylabel = "DHHN2016 Height", zlabel = "H")
end
