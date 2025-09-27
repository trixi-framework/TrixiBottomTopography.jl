# Changelog

TrixiBottomTopography.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1)
used in the Julia ecosystem. Notable changes will be documented in this file
for human readability.

## Changes in the v0.1 lifecycle

#### Added

- Implementation of visualization routines `plot_topography` and `plot_topography_with_interpolation_knots`

#### Changed

- The required Julia version is updated to v1.10.
- Visualization routines rely on Makie.jl instead of PyPlot.jl [#52]

#### Deprecated

- Helper function `evaluate_bicubicspline_interpolant` renamed to be `evaluate_two_dimensional_interpolant`

#### Removed
