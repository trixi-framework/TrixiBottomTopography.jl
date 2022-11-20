# TrixiBottomTopography.jl
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/)
[![Build Status](https://github.com/maxbertrand1996/TrixiBottomTopography.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/maxbertrand1996/TrixiBottomTopography.jl/actions/workflows/ci.yml)
[![Coverage Status](https://coveralls.io/repos/github/maxbertrand1996/TrixiBottomTopography.jl/badge.svg?branch=main)](https://coveralls.io/github/maxbertrand1996/TrixiBottomTopography.jl?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

TrixiBottomTopography.jl is a supplementary package to the numerical solver [Trixi.jl]( https://github.com/trixi-framework/Trixi.jl) which enables to use real life geographical data for the bottom topography function of the shallow water equations.

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

require a (piecewise) differentiable function $b$ which describes the bottom topography.

Geographical data is almost always given as scattered data points on a coordinate system with corresponding elevation. So to incorporate geographical data into the shallow water equations, we need to define a function, which remodels the topography from the data.

TrixiBottomTopography.jl does this by B-spline interpolation of the underlying data.

## Functionalities

This package contains the following three main functionalities:
- Converting geographical data given in form of `.xyz` files from the [DGM data set](https://www.opengeodata.nrw.de/produkte/geobasis/hm/dgm1_xyz/dgm1_xyz/) provided by [Geobasis NRW](https://www.bezreg-koeln.nrw.de/brk_internet/geobasis/hoehenmodelle/digitale_gelaendemodelle/gelaendemodell/index.html) to make it readable for `TrixiBottomTopography.jl`
- Setting up a B-spline interpolation structure in one and two dimensions which contains all the relevant information to define a B-spline interpolation function with additional specifications
- Using the B-spline structure to set up a B-spline interpolation function

A detailed description of the functionalities can be found on the following pages.

## Examples

Some examples which use the implemented functionalities, see the [examples folder](https://github.com/maxbertrand1996/TrixiBottomTopography.jl/tree/main/examples) of this repository.

## Authors
TrixiBottomTopography.jl was developed by [Maximilian Dominique Bertrand](https://github.com/maxbertrand1996) (University of Cologne, Germany).

## License and contributing
TrixiBottomTopography.jl is published under the MIT license (see [License](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/licence/)). We
are very happy to accept contributions from everyone, preferably in the form of
a PR.