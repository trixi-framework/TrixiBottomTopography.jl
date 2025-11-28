##############################################################################
# Script which uses the functionalities implemented in TrixiBottomTopography #
# to plot a one dimensional cubic B-spline interpolated section of the Rhine #
# river with free end condition.                                             #
##############################################################################

# Include packages
using TrixiBottomTopography

# Define data path
root_dir = pkgdir(TrixiBottomTopography)
data = joinpath(root_dir, "examples", "data", "rhine_data_1d_20_x.txt")

# Define B-spline structure
spline_struct = CubicBSpline(data)
# Define B-spline interpolation function
spline_func(x) = spline_interpolation(spline_struct, x)

# Evaluate the cubic B-spline on a new set of nodes and plot
if isdefined(Main, :Makie)
    # Define interpolation points
    n = 200
    x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))

    # Get interpolated values
    y_int_pts = spline_func.(x_int_pts)

    plot_topography(x_int_pts, y_int_pts;
                    xlabel = "ETRS89 East",
                    ylabel = "DHHN2016 Height")
end