# TrixiBottomTopography.jl API

```@meta
CurrentModule = TrixiBottomTopography
```

```@autodocs
Modules = [TrixiBottomTopography]
```

## Makie.jl extension

```@docs
TrixiBottomTopography.evaluate_two_dimensional_interpolant
TrixiBottomTopography.plot_topography
TrixiBottomTopography.plot_topography_with_interpolation_knots
```

## KernelInterpolation.jl extension

```@docs
TrixiBottomTopography.RBFInterpolation
TrixiBottomTopography.RBFInterpolation1D
TrixiBottomTopography.RBFInterpolation2D
```

## JuMPHiGHS.jl extension
```@docs
TrixiBottomTopography.LaverySpline1D(::Vector{T}, ::Vector{T}) where {T<:Real}
TrixiBottomTopography.LaverySpline1D(::String)
TrixiBottomTopography.LaverySpline2D(::Vector{T}, ::Vector{T}, ::Matrix{T}) where {T<:Real}
TrixiBottomTopography.LaverySpline2D(::String)
TrixiBottomTopography.spline_interpolation(::LaverySpline1D, ::Number)
TrixiBottomTopography.spline_interpolation(::LaverySpline2D, ::Number, ::Number)
```