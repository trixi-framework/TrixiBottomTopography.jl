# TrixiBottomTopography.jl

`TrixiBottomTopography.jl` is a supplementary package to the numerical solver [Trixi.jl]( https://github.com/trixi-framework/Trixi.jl), which enables to use real world geographical data for the bottom topography function of the shallow water equations.

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

require a (piecewise) differentiable function $b$, which describes the bottom topography.

Geographical data is almost always given as scattered data points on a coordinate system with the corresponding elevation. So to incorporate geographical data into the shallow water equations, we need to define a function that remodels the topography from the data.

`TrixiBottomTopography.jl` does this by B-spline interpolation of the underlying data.

## Functionalities

This package contains the following three main functionalities:
- Converting geographical data given in form of `.xyz` files from the [DGM data set](https://www.opengeodata.nrw.de/produkte/geobasis/hm/dgm1_xyz/dgm1_xyz/) provided by [Geobasis NRW](https://www.bezreg-koeln.nrw.de/brk_internet/geobasis/hoehenmodelle/digitale_gelaendemodelle/gelaendemodell/index.html) to make it readable for `TrixiBottomTopography.jl`
- Setting up a B-spline interpolation structure in one and two dimensions which contains all the relevant information to define a B-spline interpolation function with additional specifications
- Using the B-spline structure to set up a B-spline interpolation function

A detailed description of the functionalities can be found in the [documentation](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/) to this package.

## Installation
If you have not yet installed Julia, please [follow the instructions for your operating system](https://julialang.org/downloads/platform/). TrixiBottomTopography works with Julia v1.7.

As TrixiBottomTopography is **not** a registered Julia package, you have to download it locally 
```
git clone https://github.com/trixi-framework/TrixiBottomTopography.jl.git
```
and run it from within the cloned directory
```
julia --project=@.
```
to make use of the implemented functionalities.

## Examples

To see a first example of `TrixiBottomTopography.jl` a default example has been implemented. First, load the package
```julia
julia> using TrixiBottomTopography
```
Then call the `T`rixi`B`ottom`T`opograpy `default example`
```julia
julia> TBT_default_example()
```
If the implementation was successful, the following window appears:

![image](https://user-images.githubusercontent.com/101979498/203507049-279bc69b-3acc-4c55-888f-26e02c1edabe.png)

Further examples can be found in the [examples folder](https://github.com/trixi-framework/TrixiBottomTopography.jl/tree/main/examples) of this repository.

## Authors
`TrixiBottomTopography.jl` was developed by [Maximilian Dominique Bertrand](https://github.com/maxbertrand1996) (University of Cologne, Germany) with the help of [Andrew Winters](https://liu.se/en/employee/andwi94) (Linköping University, Sweden) and [Michael Schlottke-Lakemper](https://lakemper.eu/) (RWTH Aachen University, Germany).

## License and contributing
`TrixiBottomTopography.jl` is published under the MIT license (see [License](https://github.com/trixi-framework/TrixiBottomTopography.jl/blob/main/LICENSE)). We
are pleased to accept contributions from everyone, preferably in the form of
a PR.