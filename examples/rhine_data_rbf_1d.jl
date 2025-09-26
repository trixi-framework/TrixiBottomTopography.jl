###############################################################################
# Script which uses the functionalities implemented in TrixiBottomTopography  #
# to plot a one dimensional RBF interpolated section of the Rhine river with  #
# data extracted along the y-direction.                                       #
###############################################################################

# Include packages
using TrixiBottomTopography
using KernelInterpolation

# Define data path
root_dir = pkgdir(TrixiBottomTopography)
data = joinpath(root_dir, "examples", "data", "rhine_data_1d_20_x.txt")

# Define RBF structure
kernel = ThinPlateSplineKernel{1}()
itp = RBFInterpolation1D(data, kernel)
# Define interpolation function
itp_func(x) = itp(x)

# Evaluate the cubic B-spline on a new set of nodes and plot
if isdefined(Main, :Makie)
    nodes = nodeset(itp)
    x_knots = sort(unique(nodes[:, 1]))
    y_knots = itp_func.(x_knots)

    # Define interpolation points
    n = 200
    x_int_pts = Vector(LinRange(x_knots[1], x_knots[end], n))

    # Get interpolated values
    y_int_pts = itp_func.(x_int_pts)

    plot_topography_with_interpolation_knots(x_int_pts, y_int_pts, x_knots, y_knots;
                                             xlabel = "ETRS89 East",
                                             ylabel = "DHHN2016 Height",
                                             legend_position = :rb)
end
