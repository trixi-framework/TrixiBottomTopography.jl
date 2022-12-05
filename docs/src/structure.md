# B-spline interpolation structure

Once the underlying data is in the correct format, we can start defining B-spline
structures which we will use later to define the interpolation functions.

## One dimensional structures

For the one dimensional case, the structures [`LinearBSpline`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.LinearBSpline)
and [`CubicBSpline`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.CubicBSpline)
are available. They contain all relevant values to define linear and cubic B-spline
interpolation functions corresponding to linear and cubic B-spline interpolation.
These are:
- `x`: A vector of values in x-direction.
- `Delta`: The length of a single patch in the given data set. A patch is the area between two consecutive
           `x` values. `Delta` corresponds to the distance between two consecutive values in x-direction.
           As we are only considering Cartesian grids, `Delta` is equal for all patches.
- `Q`: A vector that contains the control points.
- `IP`: The coefficients matrix.

To populate the structure, the outer constructor functions [`LinearBSpline(data_path)`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.LinearBSpline-Tuple{String})
and [`CubicBSpline(data_path)`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.CubicBSpline-Tuple{String})
are implemented. These constructors use the files in `data_path` to obtain the values which
will be stored in the corresponding structure,
as well as [`LinearBSpline(x,y)`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.LinearBSpline-Tuple{Vector{T}%20where%20T,%20Vector{T}%20where%20T})
and [`CubicBSpline(x,y)`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.CubicBSpline-Tuple{Vector{T}%20where%20T,%20Vector{T}%20where%20T}) which use given vectors `x` and `y`.

To get a better idea of the constructor functions, we consider the example [rhine\_data\_cubic-nak.jl](https://github.com/trixi-framework/TrixiBottomTopography.jl/blob/main/examples/rhine_data_cubic-nak.jl)
from the [examples folder](https://github.com/trixi-framework/TrixiBottomTopography.jl/tree/main/examples) of this repo.
This particular example reads one dimensional bottom topography data from a `.txt` file
and constructs a cubic B-spline interpolation with not-a-knot end condition
and smoothing of the data.

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

For the cubic case, we can also set the optional parameters `end_condition`, which defines (as the name suggests) the end condition of the spline.
Available end conditions are the `not-a-knot` and the `free` end condition.
By default, `end_condition` is set to `free`. If you are not familiar with the differences
and influence of these end conditions, see Chapter 1 of
- Quentin Agrapart & Alain Batailly (2020), Cubic and bicubic spline interpolation in Python.
  [hal-03017566v2](https://hal.archives-ouvertes.fr/hal-03017566v2).

Besides the end condition, we can also specify a `smoothing_factor` for the cubic B-spline.
This smoothing parameter defines a trade-off for the resulting cubic B-spline interpolation
between how well the B-spline models (or fits) the original data
and minimizing the curvature by defining new `y` values.
This procedure is called **spline smoothing**. There is no general approach to determine
which `smoothing_factor` is best suited for a given problem.
It must be determined by the user via trial and error.
However, as a general rule, the larger the smoothing factor, the less curvature will be
present in the resulting B-spline.
To understand the underlying maths of the smoothing spline procedure, please see:
- Germán Rodríguez (2001),
  [Smoothing and non-parametric regression](https://docplayer.net/6006594-Smoothing-and-non-parametric-regression.html)

## Two dimensional structures

For the two dimensional case, the structures [`BilinearBSpline`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.BilinearBSpline)
and [`BicubicBSpline`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.BicubicBSpline)
are implemented, which contain all relevant values to define bilinear and bicubic B-spline interpolation functions. These are:
- `x`: Vector of values in x-direction.
- `y`: Vector of values in y-direction.
- `Delta`: Length of one side of a single patch in the given data set. A patch is the area between two consecutive `x` and `y` values. `Delta` corresponds to the distance between two consecutive values in x-direction. The implementation only considers Cartesian grids, so `Delta` is equal for all patches in x and y-direction.
- `Q`: Matrix which contains the control points.
- `IP`: Coefficients matrix.

To populate the structure, the outer constructor functions [`BilinearBSpline(data_path)`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.BilinearBSpline-Tuple{String})
and [`BicubicBSpline(data_path)`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.BicubicBSpline-Tuple{String})
are implemented. They use the files in `data_path` to obtain the values which will be stored
in the corresponding structure.
Further, [`BilinearBSpline(x,y,z)`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.BilinearBSpline-Tuple{Vector{T}%20where%20T,%20Vector{T}%20where%20T,%20Matrix{T}%20where%20T})
and [`BicubicBSpline(x,y,z)`](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.BicubicBSpline-Tuple{Vector{T}%20where%20T,%20Vector{T}%20where%20T,%20Matrix{T}%20where%20T})
which use given vectors `x` and `y` and matrix `z`.
The `x`, `y` and `z` data are organized in the following form:

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

To better understand the constructor functions, we consider the example [rhine\_data\_bicubic-nak.jl](https://github.com/trixi-framework/TrixiBottomTopography.jl/blob/main/examples/rhine_data_bicubic-nak.jl) of this repo.
This example reads in two dimensional bottom topography data from a `.txt` file and
creates a bicubic B-spline interpolation with not-a-knot end condition
and smoothing of the data.

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

For the bicubic spline, we can set the optional parameters `end_condition`, which defines (as in the one dimensional case) the end condition of the spline.
Again, the available end conditions are the `not-a-knot` and the `free` end condition.
By default, `end_condition` is set to `free`.
If you are not familiar with the differences and influence of these end conditions,
see Chapter 2 of
- Quentin Agrapart & Alain Batailly (2020), Cubic and bicubic spline interpolation in Python.
  [hal-03017566v2](https://hal.archives-ouvertes.fr/hal-03017566v2)

Besides the end condition, we can also specify a `smoothing_factor` for the bicubic spline.
Just as in the one dimensional constructor, this smoothing provides a trade-off
of the bicubic B-spline interpolation between how well it models (or fits) the
original topography data and minimizing the curvature of the spline
by defining new `z` values.
This procedure is called **thin plate spline**. There is no general approach to which `smoothing_factor` is best suited for the problem, and it must be determined via
trial and error by the user. To understand the underlying maths of the thin plate spline,
please see:
- Gianluca Donato and Serge Belongie (2001),
  Approximate Thin Plate Spline Mappings
  [DOI: 10.1007/3-540-47977-5_2](https://link.springer.com/content/pdf/10.1007/3-540-47977-5_2.pdf)