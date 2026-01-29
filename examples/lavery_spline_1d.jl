##############################################################################
# Script which uses the functionalities implemented in TrixiBottomTopography #
# to plot a one dimensional Lavery spline of data with sharp transitions.    #
##############################################################################

# Include packages
using TrixiBottomTopography

# Load in data from one of the examples from Eriksson and Jemsson thesis
# https://liu.diva-portal.org/smash/get/diva2:1918338/FULLTEXT01.pdf
xData = Vector(0:29)
yData = [0, 0, 0, 0, 0, 4, 2, 4, 3, 2,
    1, 0, 9, 9, 9, 0, 0, 0, 6, 7,
    5, 6, 7, 8, 0, 0, 0, 0, 0, 0]

# Lavery spline that does not generate new extrema
spline_struct = LaverySpline1D(xData, yData)
# Define Lavery spline interpolation function
spline_func(x) = spline_interpolation(spline_struct, x)

# Evaluate the cubic Lavery spline on a new set of nodes and plot
if isdefined(Main, :Makie)
    # Define interpolation points
    n = 500
    x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))

    # Get interpolated values
    y_int_pts = spline_func.(x_int_pts)

    # Get the original interpolation knots
    x_knots = spline_struct.x
    y_knots = spline_func.(x_knots)

    plot_topography_with_interpolation_knots(x_int_pts, y_int_pts, x_knots, y_knots;
                                             xlabel = "x",
                                             ylabel = "y",
                                             legend_position = :lt)
end
