##############################################################################
# Script which uses the functionalities implemented in TrixiBottomTopography #
# to plot a bicubic B-spline interpolated section of the Rhine river with    #
# not-a-knot end condition.                                                  #
##############################################################################

# Include packages
using TrixiBottomTopography
using Plots

# Helperfunction to fill the solution matrix
# Input parameters:
#  - f: spline function
#  - x: vector of x values
#  - y: vector of y values
function fill_sol_mat(f, x, y)

  # Get dimensions for solution matrix
  n = length(x)
  m = length(y)

  # Create empty solution matrix
  z = zeros(n,m)

  # Fill solution matrix
  for i in 1:n, j in 1:m
    # Evaluate spline functions
    # at given x,y values
    z[j,i] = f(x[i], y[j])
  end

  # Return solution matrix
  return z
end

# Define data path
data = joinpath(@__DIR__, "data", "rhine_data_2d_20.txt")

# Define B-spline structure
spline_struct = BicubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 9999)
# Define B-spline interpolation function
spline_func(x,y) = spline_interpolation(spline_struct, x, y)

# Define interpolation points
n = 100
x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))
y_int_pts = Vector(LinRange(spline_struct.y[1], spline_struct.y[end], n))

# Get interpolated matrix
z_int_pts = fill_sol_mat(spline_func, x_int_pts, y_int_pts)

# Plotting
pyplot()
surface(x_int_pts, y_int_pts, z_int_pts, camera=(-30,30),
        xlabel="E", ylabel="N", zlabel="H",
        label="Bottom topography",
        title="Bicubic B-spline interpolation\nwith not-a-knot end condition and smoothing")