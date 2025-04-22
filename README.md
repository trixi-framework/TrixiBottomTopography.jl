# TrixiBottomTopography.jl
[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://trixi-framework.github.io/TrixiBottomTopography.jl/stable)
[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/)
[![Build Status](https://github.com/trixi-framework/TrixiBottomTopography.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/trixi-framework/TrixiBottomTopography.jl/actions/workflows/ci.yml)
[![Coverage Status](https://coveralls.io/repos/github/trixi-framework/TrixiBottomTopography.jl/badge.svg?branch=main)](https://coveralls.io/github/trixi-framework/TrixiBottomTopography.jl?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15122147.svg)](https://doi.org/10.5281/zenodo.15122147)

`TrixiBottomTopography.jl` is a supplementary package to the numerical solvers [Trixi.jl](https://github.com/trixi-framework/Trixi.jl) and [TrixiShallowWater.jl](https://github.com/trixi-framework/TrixiShallowWater.jl), which enables use of real world geographical data for the bottom topography function of the shallow water equations.

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
- Converting geographical data given in form of `.xyz` files, e.g., from the [DGM data set](https://www.opengeodata.nrw.de/produkte/geobasis/hm/) provided by [Geobasis NRW](https://www.bezreg-koeln.nrw.de/geobasis-nrw) to make it readable for `TrixiBottomTopography.jl`
- Setting up a B-spline interpolation structure in one and two dimensions which contains all the relevant information to define a B-spline interpolation function with additional specifications
- Using the B-spline structure to set up a B-spline interpolation function

A detailed description of the functionalities can be found in the [documentation](https://trixi-framework.github.io/TrixiBottomTopography.jl/dev/) to this package.

## Installation
If you have not yet installed Julia, please [follow the instructions for your
operating system](https://julialang.org/downloads/platform/). TrixiBottomTopography.jl works
with Julia v1.10 and newer. We recommend using the latest stable release of Julia.

### For users
TrixiBottomTopography.jl and its related tools are registered Julia packages. Hence, you
can install it by executing the following commands in the Julia REPL:
```julia
julia> import Pkg; Pkg.add("TrixiBottomTopography")
```

The available visualization functionality uses [Makie.jl](https://github.com/JuliaPlots/Makie.jl/).
A Makie backend, such as [GLMakie.jl](https://github.com/JuliaPlots/GLMakie.jl/), can
be installed in addition to TrixiBottomTopography
```julia
julia> using Pkg; Pkg.add("GLMakie")
```

To use TrixiBottomTopography.jl together with the numerical solver framework [Trixi.jl](https://github.com/trixi-framework/Trixi.jl),
you need both Trixi.jl and a relevant time integration sub-package of
[OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl), e.g.,
for high-order low-storage Runge-Kutta schemes. These can be added
by executing
```julia
julia> using Pkg; Pkg.add(["Trixi", "OrdinaryDiffEqLowStorageRK"])
```
TrixiBottomTopography.jl can also be used together with
[TrixiShallowWater.jl](https://github.com/trixi-framework/TrixiShallowWater.jl), a solver suite specifically designed for shallow water flow applications.
An example that combines TrixiBottomTopography.jl with wet/dry transitions and
shock capturing to model a tsunami runup is available as a
[tutorial](https://trixi-framework.github.io/TrixiShallowWater.jl/stable/tutorials/elixir_shallowwater_monai_tsunami/)
in TrixiShallowWater.jl.

### For developers
If you plan on editing TrixiBottomTopography.jl itself, you can download TrixiBottomTopography.jl
locally and use the code from the cloned directory:
```bash
git clone git@github.com:trixi-framework/TrixiBottomTopography.jl.git
cd TrixiBottomTopography.jl
mkdir run
cd run
julia --project=. -e 'using Pkg; Pkg.develop(PackageSpec(path=".."))' # Install local TrixiBottomTopography.jl clone
julia --project=. -e 'using Pkg; Pkg.add(["GLMakie", "Trixi", "OrdinaryDiffEqLowStorageRK"])' # Install additional packages
```
Note that the additional packages are optional and can be omitted.

If you installed TrixiBottomTopography.jl this way, you always have to start Julia with the `--project`
flag set to your `run` directory, e.g.,
```bash
julia --project=.
```
if already inside the `run` directory or
```bash
julia --project=run
```
if inside the repository root directory.

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
The full list of contributors can be found in [AUTHORS.md](AUTHORS.md).

## License and contributing
`TrixiBottomTopography.jl` is published under the MIT license (see [LICENSE.md](LICENSE.md)).
We are pleased to accept contributions from everyone, preferably in the form of a PR.