# B-spline interpolation function

If the B-spline structure is defined, it can be used to define the B-spline interpolation function. In this chapter, we are going to continue with the examples from the previous chapter for the one and two dimensional case.

## One dimensional case

In the previous chapter we started looking at the example file [rhine\_data\_cubic-nak.jl](https://github.com/maxbertrand1996/TrixiBottomTopography.jl/blob/9f6c7e967a3b094dbfa43688d25a8998fce40014/examples/rhine_data_cubic-nak.jl) from the [examples folder](https://github.com/maxbertrand1996/TrixiBottomTopography.jl/tree/9f6c7e967a3b094dbfa43688d25a8998fce40014/examples) where we already defined the B-spline structure for the cubic B-spline interpolation with not-a-knot end condition and smoothing.

```julia
# Include packages
using TrixiBottomTopography
using Plots

# Get root directory
dir_path = pkgdir(TrixiBottomTopography)

# Define data path
data = string(dir_path, "/examples/data/rhine_data_1d_20_x.txt")

# Define B-spline structure
spline_struct = CubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 999)
```

To define the B-spline interpolation function with respect to a variable `x`, we use the `spline_interpolation` function which for the one dimensional case which is implemented as [`spline_interpolation(linear_struct, x)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.spline_interpolation-Tuple{LinearBSpline,%20Any}) and [`spline_interpolation(cubic_struct, x)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.spline_interpolation-Tuple{CubicBSpline,%20Any}).

```julia
# Define B-spline interpolation function
spline_func(x) = spline_interpolation(spline_struct, x)
```

This defines the cubic B-spline interpolation function with not-a-knot end condition and smoothing with respect to variable `x` because of how we defined `spline_struct` previously. If we want to visualize the interpolation function with 100 interpolation points, we define:
```julia
# Define interpolation points
n = 100
x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))
```

and receive the corresponding `y` values by:

```julia
# Get interpolated values
y_int_pts = spline_func.(x_int_pts)
```

Plotting the interpolated points, using 

```julia
# Plotting
pyplot()
plot(x_int_pts, y_int_pts,
     xlabel="ETRS89 East", ylabel="DHHN2016 Height",
     label="Bottom topography", 
     title="Cubic B-spline interpolation with not-a-knot end condition and smoothing")
```

gives the following representation:

![image](https://raw.githubusercontent.com/maxbertrand1996/TrixiBottomTopography.jl/main/assets/images/Cubic_nak_smth.png)

## Two dimensional case

In the previous chapter we started looking at the example file [rhine\_data\_bicubic-nak.jl](https://github.com/maxbertrand1996/TrixiBottomTopography.jl/blob/main/examples/rhine_data_bicubic-nak.jl) from the [examples folder](https://github.com/maxbertrand1996/TrixiBottomTopography.jl/tree/9f6c7e967a3b094dbfa43688d25a8998fce40014/examples) where we already defined the B-spline structure for the bicubic B-spline interpolation with not-a-knot end condition and smoothing.

```julia
# Include packages
using TrixiBottomTopography
using Plots

# Get root directory
dir_path = pkgdir(TrixiBottomTopography)

# Define data path
data = string(dir_path, "/examples/data/rhine_data_2d_20.txt")

# Define B-spline structure
spline_struct = BicubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 9999)
```

To define the B-spline interpolation function with respect to a variables `x` and `y`, we use the `spline_interpolation` function which for the two dimensional case which is implemented as [`spline_interpolation(bilinear_struct, x, y)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.spline_interpolation-Tuple{BilinearBSpline,%20Any,%20Any}) and [`spline_interpolation(bicubic_struct, x, y)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.spline_interpolation-Tuple{BicubicBSpline,%20Any,%20Any}).

```julia
# Define B-spline interpolation function
spline_func(x,y) = spline_interpolation(spline_struct, x, y)
```
This defines the bicubic B-spline interpolation function with not-a-knot end condition and smoothing with respect to variables `x` and `y` because of how we defined `spline_struct` previously. If we want to visualize the interpolation function with 100 interpolation points in each direction, we define:

```julia
# Define interpolation points
n = 100
x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))
y_int_pts = Vector(LinRange(spline_struct.y[1], spline_struct.y[end], n))
```

To fill a matrix `z_int_pts` which contains the corresponding `z` values for `x_int_pts` and `y_int_pts`, we define a helperfunction `fill_sol_mat`:

```julia
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
```

and receive the `z_int_pts` values by setting:

```julia
# Get interpolated matrix
z_int_pts = fill_sol_mat(spline_func, x_int_pts, y_int_pts)
```

Plotting the interpolated values

```julia
# Plotting
pyplot()
surface(x_int_pts, y_int_pts, z_int_pts, camera=(-30,30),
        xlabel="ETRS89\n East", ylabel="ETRS89\n North", zlabel="DHHN2016\n Height", 
        label="Bottom topography", 
        title="Cubic B-spline interpolation with not-a-knot end condition")
```

gives the following representation:

![image](https://raw.githubusercontent.com/maxbertrand1996/TrixiBottomTopography.jl/main/assets/images/Bicubic_nak_smth.png)