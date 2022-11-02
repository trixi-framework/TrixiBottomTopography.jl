##############################################################################
# Script which uses the functionalities implemented in TrixiBottomTopography #
# to polt a bicubic B-spline interpolated section of the Rhine river with    #
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

# Get root directory
dir_path = pkgdir(TrixiBottomTopography)

# Define data path
data = string(dir_path, "/examples/data/rhine_data_2d_20.txt")

# Define B-spline structure
spline_struct = BicubicBSpline(data; end_condition = "not-a-knot")
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
        xlabel="ETRS89\n East", ylabel="ETRS89\n North", zlabel="DHHN2016\n Height", 
        label="Bottom topography", 
        title="Cubic B-spline interpolation with not-a-knot end condition")