# B-spline interpolation function

If the B-spline structure is defined, it can be used to define the B-spline interpolation function.
In this chapter, we will continue with the examples from the [B-spline structure](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/structure/)
section for the one and two dimensional case.

## One dimensional case

In the [B-spline structure](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/structure/)
section, we began with the example file
[rhine\_data\_cubic-nak.jl](https://github.com/trixi-framework/TrixiBottomTopography.jl/blob/main/examples/rhine_data_cubic-nak.jl)
from the [examples folder](https://github.com/trixi-framework/TrixiBottomTopography.jl/tree/main/examples)
where we already defined the B-spline structure for the cubic B-spline interpolation with not-a-knot end condition and smoothing.

```julia
# Include packages
using TrixiBottomTopography
using CairoMakie

# Define data path
root_dir = pkgdir(TrixiBottomTopography)
data = joinpath(root_dir, "examples", "data", "rhine_data_1d_20_x.txt")

# Define B-spline structure
spline_struct = CubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 999)
```

To define the B-spline interpolation function for a variable `x`,
we use the `spline_interpolation` function. For the one dimensional case,
this function is implemented as [`spline_interpolation(linear_struct, x)`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.spline_interpolation-Tuple{LinearBSpline,%20Any})
and [`spline_interpolation(cubic_struct, x)`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.spline_interpolation-Tuple{CubicBSpline,%20Any}).

```julia
# Define B-spline interpolation function
spline_func(x) = spline_interpolation(spline_struct, x)
```

This defines the cubic B-spline interpolation function with not-a-knot end condition
and smoothing with respect to variable `x` from the previously created `spline_struct`.
If we want to visualize the interpolation function with 100 interpolation points, we define the following:
```julia
# Define interpolation points
n = 200
x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))
```

and evaluate to obtain the corresponding `y` values by:

```julia
# Get interpolated values
y_int_pts = spline_func.(x_int_pts)
```

Plotting the interpolated points can be done via

```julia
plot_topography(x_int_pts, y_int_pts; xlabel = "ETRS89 East", ylabel = "DHHN2016 Height")
```

which gives the following representation:

![image](https://github.com/user-attachments/assets/ceb242a8-53a6-424f-abc9-39789d6e31d4)

Alternatively, one can plot the interpolated bottom topography together
with the interpolation knots as follows

```julia
# Get the original interpolation knots
x_knots = spline_struct.x
y_knots = spline_func.(x_knots)

plot_topography_with_interpolation_knots(x_int_pts, y_int_pts, x_knots, y_knots;
                                         xlabel = "ETRS89 East", ylabel = "DHHN2016 Height" )
```

which gives

![image](https://github.com/user-attachments/assets/91c515f8-986b-42ff-bf5d-1927687a59c5)

## Two dimensional case

In the [B-spline structure](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/structure/)
section, we examined the example file
[rhine\_data\_bicubic-nak.jl](https://github.com/trixi-framework/TrixiBottomTopography.jl/blob/main/examples/rhine_data_bicubic-nak.jl)
from the [examples folder](https://github.com/trixi-framework/TrixiBottomTopography.jl/tree/main/examples)
where we already created the B-spline structure for the bicubic B-spline interpolation
with not-a-knot end condition and smoothing.

```julia
# Include packages
using TrixiBottomTopography
using CairoMakie

# Define data path
root_dir = pkgdir(TrixiBottomTopography)
data = joinpath(root_dir, "examples", "data", "rhine_data_2d_20.txt")

# Define B-spline structure
spline_struct = BicubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 9999)
```

To define the B-spline interpolation function for variables `x` and `y`,
we use the `spline_interpolation` function for the two dimensional case.
This functionality is implemented as [`spline_interpolation(bilinear_struct, x, y)`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.spline_interpolation-Tuple{BilinearBSpline,%20Any,%20Any})
and [`spline_interpolation(bicubic_struct, x, y)`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.spline_interpolation-Tuple{BicubicBSpline,%20Any,%20Any}).

```julia
# Define B-spline interpolation function
spline_func(x,y) = spline_interpolation(spline_struct, x, y)
```
This defines the bicubic B-spline interpolation function with not-a-knot end condition
and smoothing for variables `x` and `y` due to the previously constructed `spline_struct`. If we want to visualize the bicubic interpolation function with 100
interpolation points in each spatial direction, we define:

```julia
# Define interpolation points
n = 100
x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))
y_int_pts = Vector(LinRange(spline_struct.y[1], spline_struct.y[end], n))
```

To fill a matrix `z_int_pts`, which contains the corresponding `z` values
for `x_int_pts` and `y_int_pts`, we use the helper function
`evaluate_bicubicspline_interpolant` implemented in `ext/TrixiBottomTopographyMakieExt.jl`.

```julia
# Get interpolated matrix
z_int_pts = evaluate_bicubicspline_interpolant(spline_func, x_int_pts, y_int_pts)
```

Plotting the interpolated values

```julia
  plot_topography(x_int_pts, y_int_pts, z_int_pts;
                  xlabel="ETRS89\n East",
                  ylabel="ETRS89\n North",
                  zlabel="DHHN2016\n Height",
                  azimuth_angle = 54 * pi / 180,
                  elevation_angle = 27 * pi / 180)
```

gives the following representation:
TODO: remake figure
![image](https://github.com/user-attachments/assets/1203483a-b414-45b1-a69a-c6e284eeb0c2)

Alternatively, one can plot the interpolated bottom topography together
with the interpolation knots as follows

```julia
# Get the original interpolation knots
x_knots = spline_struct.x
y_knots = spline_struct.y
z_knots = evaluate_bicubicspline_interpolant(spline_func, x_knots, y_knots)

plot_topography_with_interpolation_knots(x_int_pts, y_int_pts, z_int_pts,
                                         x_knots, y_knots, z_knots;
                                         xlabel="ETRS89\n East",
                                         ylabel="ETRS89\n North",
                                         zlabel="DHHN2016\n Height",
                                         azimuth_angle = 54 * pi / 180,
                                         elevation_angle = 27 * pi / 180)
```

which gives
![image](https://github.com/user-attachments/assets/ed082318-47ec-4680-8cd0-04997c3a95b4)