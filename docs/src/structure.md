# B-spline interpolation structure
If you have the underlying data in the correct format, you can start defining B-spline structures which are later used to define the interpolation functions. 

## One dimensional structures

For the one dimensional case, the structures [`LinearBSpline`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.LinearBSpline) and [`CubicBspline`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.CubicBSpline) are implemented which contain all relevant values to define linear and cubic B-spline interpolation functions correspond to linear and cubic B-spline interpolation. These are:
- `x`: A vector of values in x-direction
- `h`: The length of a single patch in the given data set. A patch is the area between two consecutive 
       `x` values. `h` corresponds to the distance between two consecutive values in x-direction. 
       As we are only considering Cartesian grids, `h` is equal for all patches
- `Q`: A vector which contains the control points
- `IP`: The coefficients matrix

To populate the structure, the outer constructor functions [`LinearBSpline(data_path)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.LinearBSpline-Tuple{String}) and [`CubicBSpline(data_path)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.CubicBSpline-Tuple{String}) which use the files in `data_path` to obtain the values which will be stored in the corresponding structure, as well as [`LinearBSpline(x,y)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.LinearBSpline-Tuple{Vector{T}%20where%20T,%20Vector{T}%20where%20T}) and [`CubicBSpline(x,y)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.CubicBSpline-Tuple{Vector{T}%20where%20T,%20Vector{T}%20where%20T}) which use given vectors `x` and `y`.

To get a better idea of the constructor functions, we are taking a look at example [rhine\_data\_cubic-nak.jl](https://github.com/maxbertrand1996/TrixiBottomTopography.jl/blob/9f6c7e967a3b094dbfa43688d25a8998fce40014/examples/rhine_data_cubic-nak.jl) from the [examples folder](https://github.com/maxbertrand1996/TrixiBottomTopography.jl/tree/9f6c7e967a3b094dbfa43688d25a8998fce40014/examples) of this repo which reads one dimensional bottom topography data from a `.txt` file and does a cubic B-spline interpolation with not-a-knot end condition and smoothing of the data.

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

For the cubic case we can also set the optional parameters `end_condition` which defines (as the name suggests) the end condition of the spline. Implemented are the `not-a-knot` and the `free` end condition. By default, `end_condition` is set to `free`. If you are not familiar with the differences between these end conditions, see Chapter 1 of
- Quentin Agrapart & Alain Batailly (2020),
Cubic and bicubic spline interpolation in Python. 
[hal-03017566v2](https://hal.archives-ouvertes.fr/hal-03017566v2)

Besides the end condition, in the cubic case, we can also specify a `smoothing_factor`, which defines trade-off degree of the cubic B-spline interpolation of goodness of fit and minimizing the curvature by defining new `y` values. This proceedure is called **spline smooothing**. There is no general approach which `smoothing_factor` is best suited for the problem and has to be determined via trial and error. To understand the underlying maths, please see:
- Germán Rodríguez (2001),
  [Smoothing and non-parametric regression](https://data.princeton.edu/eco572/smoothing.pdf)

In the next chapter, we are going to continue with the presented example.

## Two dimensional structures

For the one dimensional case, the structures [`BilinearBSpline`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.BilinearBSpline) and [`BicubicBSpline`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.BicubicBSpline) are implemented which contain all relevant values to define bilinear and bicubic B-spline interpolation functions correspond to linear and cubic B-spline interpolation. These are:

- `x`: Vector of values in x-direction
- `y`: Vector of values in y-direction
- `h`: Length of one side of a single patch in the given data set. A patch is the area between two 
       consecutive `x` and `y` values. `h` corresponds to the distance between two consecutive 
       values in x-direction. As we are only considering Cartesian grids, `h` is equal for all 
       patches in x and y-direction
- `Q`: Matrix which contains the control points
- `IP`: Coefficients matrix

To populate the structure, the outer constructor functions [`BilinearBSpline(data_path)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.BilinearBSpline-Tuple{String}) and [`BicubicBSpline(data_path)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.BicubicBSpline-Tuple{String}) which use the files in `data_path` to obtain the values which will be stored in the corresponding structure, as well as [`BilinearBSpline(x,y,z)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.BilinearBSpline-Tuple{Vector{T}%20where%20T,%20Vector{T}%20where%20T,%20Matrix{T}%20where%20T}) and [`BicubicBSpline(x,y,z)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.BicubicBSpline-Tuple{Vector{T}%20where%20T,%20Vector{T}%20where%20T,%20Matrix{T}%20where%20T}) which use given vectors `x` and `y` and the matrix `z` which has the following form:

```math
\begin{aligned}
\begin{matrix}
    & & x_1 & x_2 & ... & x_n\\
    & & & & &\\
    y_1 & & z_{11} & z_{12} & ... & z_{1n}\\
    y_1 & & z_{21} & z_{22} & ... & z_{2n}\\
    \vdots & & \vdots & \vdots & \ddots & \vdots\\
    y_m & & z_{m1} & z_{m2} & ... & z_{mn}
  \end{matrix}
\end{aligned}
```

To get a better idea of the constructor functions, we are taking a look at example [rhine\_data\_bicubic-nak.jl](https://github.com/maxbertrand1996/TrixiBottomTopography.jl/blob/main/examples/rhine_data_bicubic-nak.jl) of this repo which reads two dimensional bottom topography data from a `.txt` file and does a bicubic B-spline interpolation with not-a-knot end condition and smoothing of the data.

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

For the bicubic case we can also set the optional parameters `end_condition` which defines (as the name suggests) the end condition of the spline. As in the one dimensional case, there are  the `not-a-knot` and the `free` end condition implemented. By default, `end_condition` is set to `free`. If you are not familiar with the differences between these end conditions, see Chapter 2 of
- Quentin Agrapart & Alain Batailly (2020),
Cubic and bicubic spline interpolation in Python. 
[hal-03017566v2](https://hal.archives-ouvertes.fr/hal-03017566v2)

Besides the end condition, in the bicubic case, we can also specify a `smoothing_factor`, which defines trade-off degree of the cubic B-spline interpolation of goodness of fit and minimizing the curvature by defining new `y` values. This procedure is called **thin plate spline**. There is no general approach which `smoothing_factor` is best suited for the problem and has to be determined via trial and error. To understand the underlying maths, please see:
- Gianluca Donato and Serge Belongie (2001),
  Approximate Thin Plate Spline Mappings
  [DOI: 10.1007/3-540-47977-5_2](https://link.springer.com/content/pdf/10.1007/3-540-47977-5_2.pdf)

In the next chapter, we are going to continue with the presented example.