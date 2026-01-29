##############################################################################
# Script which uses the functionalities implemented in TrixiBottomTopography #
# to plot a two dimensional Lavery spline of data with sharp transitions.    #
##############################################################################

# Include packages
using TrixiBottomTopography

# Load in data from one of the examples from Eriksson and Jemsson thesis
# https://liu.diva-portal.org/smash/get/diva2:1918338/FULLTEXT01.pdf
x_knots = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
y_knots = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
z_knots = [1 1 1 0 0 4 4 4 4 1;
           1 1 1 0 0 1 1 1 1 1;
           1 1 1 3 3 3 3 3 3 3;
           2 2 2 1 0 0 0 0 0 0;
           3 3 3 1 0 0 0 0 0 0;
           4 4 4 1 0 0 2 2 0 0;
           1 0 3 1 0 0 2 2 0 0;
           1 0 3 1 0 0 0 0 0 0;
           1 0 3 0 0 0 0 0 0 0;
           0 0 0 0 0 0 0 0 0 0]

# Lavery spline that does not generate new extrema
spline_struct = LaverySpline2D(x_knots, y_knots, z_knots)
# Define Lavery spline interpolation function
spline_func(x, y) = spline_interpolation(spline_struct, x, y)

# Evaluate the bicubic Lavery spline on a new set of nodes and plot
if isdefined(Main, :Makie)
    # Define interpolation points
    n = 150
    x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))
    y_int_pts = Vector(LinRange(spline_struct.y[1], spline_struct.y[end], n))

    # Get interpolated matrix
    z_int_pts = evaluate_two_dimensional_interpolant(spline_func, x_int_pts, y_int_pts)

    plot_topography_with_interpolation_knots(x_int_pts, y_int_pts, z_int_pts,
                                             x_knots, y_knots, z_knots;
                                             xlabel = "x",
                                             ylabel = "y",
                                             zlabel = "z")
end
