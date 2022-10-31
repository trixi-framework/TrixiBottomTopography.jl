# TrixiBottomTopography.jl
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/)
[![Build Status](https://github.com/maxbertrand1996/TrixiBottomTopography.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/maxbertrand1996/TrixiBottomTopography.jl/actions/workflows/ci.yml)
[![Coverage Status](https://coveralls.io/repos/github/maxbertrand1996/TrixiBottomTopography.jl/badge.svg?branch=main)](https://coveralls.io/github/maxbertrand1996/TrixiBottomTopography.jl?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

TrixiBottomTopography.jl is a supplementary package to the numerical solver [Trixi.jl]( https://github.com/trixi-framework/Trixi.jl) which enables to use real life geographical data for the bottom topography function of the shallow water equations.

## Introduction
The shallow water equations in one dimension
```math
\begin{aligned}
\begin{pmatrix} h \\ hv \end{pmatrix}_t 
+ \begin{pmatrix} hv \\ hv^2 + \frac{1}{2}gh^2 \end{pmatrix}_x 
= \begin{pmatrix} 0\\ -ghb_x \end{pmatrix}
\end{aligned}
```
and two dimensions
```math
\begin{aligned}
\begin{pmatrix} h \\ hv_1 \\ hv_2 \end{pmatrix}_t 
+ \begin{pmatrix} hv_1 \\ hv_1^2 + \frac{1}{2}gh^2 \\ hv_1v_2 \end{pmatrix}_x
+ \begin{pmatrix} hv_2 \\ hv_1v_2 \\ hv_2^2 + \frac{1}{2}gh^2  \end{pmatrix}_y
= \begin{pmatrix} 0\\ -ghb_x \\ -ghb_y \end{pmatrix}
\end{aligned}
```

require a differentiable function $b$ which describes the bottom topography.

Geographical data is almost always given as scattered data points on a coordinate system with corresponding elevation. So to incorporate geographical data into the shallow water equations, we need to define a function, which remodels the topography from the data. 

TrixiBottomTopography.jl does this by B-spline interpolation of the underlying data.

## Overview
The functionalities of this package can be summarized into the following three:
- Converting geographical data from the [DGM NRW data set](https://www.opengeodata.nrw.de/produkte/geobasis/hm/dgm1_xyz/dgm1_xyz/)
- Defining a B-spline interpolation structure
- Defining the B-spline interpolation function
### Converting geographical data
[Geobasis NRW](https://www.bezreg-koeln.nrw.de/brk_internet/geobasis/hoehenmodelle/digitale_gelaendemodelle/gelaendemodell/index.html) provides a [geographical data set](https://www.opengeodata.nrw.de/produkte/geobasis/hm/dgm1_xyz/dgm1_xyz/) of the whole German state of Nort Rhine-Westphalia. This data set contains of patches of $1km^2$ where each patch has the elevation data for 1.000.000 data points equally distributed as a grid with grid size of $1m$. The data is organized as follows:
```
357000.00 5646999.00 47.40 
357001.00 5646999.00 47.43 
357002.00 5646999.00 47.49 
357003.00 5646999.00 47.47 
357004.00 5646999.00 47.39 
357005.00 5646999.00 47.30 
357006.00 5646999.00 47.24 
...       ...        ...
```
where the first column provides the corresponding ETRS89 East coordinates, the second column the ETRS89 North coordinates and the third column the DHHN2016 height.

This data format is not accepted by the functionalities implemented in TrixiBottomTopography.jl. To make the data set readable for this package, the functions `convert_dgm_1d` and `convert_dgm_2d` are provided.

The one dimensional data has to be of the following form:
```
# Number of x values
n
# x values
x_1
...
x_n
# y values
y_1
...
y_n
```
Executing
```julia
julia> convert_dgm_1d(path_read, path_write)
```
converts the DGM data in `path_read` to the accepted data format and saves it into `path_write`. To get more information on this function, see [`convert_dgm_1d`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.convert_dgm_1d-Tuple{String,%20String}).

The two dimensional data has to be of the following form:
```
# Number of x values
n
# Number of y values
m
# x values
x_1
...
x_n
# y values
y_1
...
y_m
# z values
z_1,1
z_1,2
...
z_1,n
z_2,1
...
z_m,n
```
Executing
```julia
julia> convert_dgm_2d(path_read, path_write)
```
converts the DGM data in `path_read` to the accepted data format and saves it into `path_write`. To get more information on this function, see [`convert_dgm_2d`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.convert_dgm_2d-Tuple{String,%20String}).

### B-spline interpolation structure
If you have the underlying data in the correct format, you can start defining B-spline structures which are later used to define the interpolation functions. 

For the one dimensional case, the structures [`LinearBSpline`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.LinearBSpline) and [`CubicBspline`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.CubicBSpline) are implemented. Which correspond to linear and cubic B-spline interpolation. To populate linear structure, call the following function:
```julia
julia> spline_struct = linear_b_spline(data)
```
Assume that `data` contains a path to a data file with the values you want to interpolate in the correct format. It is also possible to directly pass vectors of `x` and `y` values to the function:
```julia
julia> spline_struct = linear_b_spline(x,y)
```
For more information see [`linear_b_spline(data)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.linear_b_spline-Tuple{String}) and [`linear_b_spline(x,y)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.linear_b_spline-Tuple{Vector{T}%20where%20T,%20Vector{T}%20where%20T}).

Populating the cubic structure is analogous to the linear case, using the function `cubic_b_spline`. For the cubic case we can additionally specify
-  The end condition
   - Free end condition
   - Not-a-knot end condition
- A smoothing factor

For further information, see [`cubic_b_spline(data)`]() and [`cubic_b_spline(x,y)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.cubic_b_spline-Tuple{Vector{T}%20where%20T,%20Vector{T}%20where%20T})

For the two dimensional case, the structures [`BilinearBSpline`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.BilinearBSpline) and [`BicubicBSpline`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.BicubicBSpline) are implemented.

Populating these works analogously to the one dimensional case with the functions
- [`bilinear_b_spline(data)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.bilinear_b_spline-Tuple{String})
- [`bilinear_b_spline(x,y,z)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.bilinear_b_spline-Tuple{Vector{T}%20where%20T,%20Vector{T}%20where%20T,%20Matrix{T}%20where%20T})
- [`bicubic_b_spline(data)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.bicubic_b_spline-Tuple{String})
- [`bicubic_b_spline(x,y,z)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.bicubic_b_spline-Tuple{Vector{T}%20where%20T,%20Vector{T}%20where%20T,%20Matrix{T}%20where%20T})

Note that `z` has to be a matrix of the following form

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


### B-spline interpolation function

With the defined structures, we can use these to set up the B-spline interpolation functions. The function to do so uses function overloading and therefore has the same name for all cases. The cases are distinguished by the given structure.

1D:
```julia
julia> spline_func(x) = spline_interpolation(spline_struct, x)
```
2D:
```julia
julia> spline_func(x, y) = spline_interpolation(spline_struct, x, y)
```

For further information, see:
- [`spline_interpolation(linear)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.spline_interpolation-Tuple{TrixiBottomTopography.LinearBSpline,%20Any}) 
- [`spline_interpolation(cubic)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.spline_interpolation-Tuple{TrixiBottomTopography.CubicBSpline,%20Any})
- [`spline_interpolation(bilinear)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.spline_interpolation-Tuple{TrixiBottomTopography.BilinearBSpline,%20Any,%20Any})
- [`spline_interpolation(bicubic)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.spline_interpolation-Tuple{TrixiBottomTopography.BicubicBSpline,%20Any,%20Any})

## Examples

For some further examples on how to use the functionalities, see the [examples folder](https://github.com/maxbertrand1996/TrixiBottomTopography.jl/tree/main/examples) of this repository.

## Authors
TrixiBottomTopography.jl was developed and is maintained by [Maximilian Dominique Bertrand](https://github.com/maxbertrand1996) (University of Cologne, Germany). 

## License and contributing
TrixiBottomTopography.jl is published under the MIT license (see [License](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/licence/)). We
are very happy to accept contributions from everyone, preferably in the form of
a PR.
