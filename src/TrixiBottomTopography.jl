"""
    TrixiBottomTopography

**TrixiBottomTopography.jl** is a supporting framework for Trixi.jl and
TrixiShallowWater.jl, which can be used to approximate bottom topography
functions using B-splines from real world data.
"""
module TrixiBottomTopography

# Include necessary packages
using LinearAlgebra: norm, diagm, qr, Tridiagonal, SymTridiagonal
using SparseArrays: sparse, spzeros
using StaticArrays: SVector, @SVector, SMatrix, @SMatrix

# Include one dimensional B-spline interpolation
include("1D/spline_utils_1D.jl")
include("1D/spline_cache_1D.jl")
include("1D/spline_methods_1D.jl")

# Include two dimensional B-spline interpolation
include("2D/spline_utils_2D.jl")
include("2D/spline_cache_2D.jl")
include("2D/spline_methods_2D.jl")

# Include auxiliary functions
include("auxiliary/convert.jl")
include("auxiliary/default_example.jl")

# Export the functions which are used for B-spline interpolation
export LinearBSpline, CubicBSpline
export BilinearBSpline, BicubicBSpline
export spline_interpolation

# Export the functions which are used DGM data conversion
export convert_dgm_1d, convert_dgm_2d

# Note, empty routines for visualization are included and exported. They are extended
# in `ext/MakieExt.jl` where their implementations are found.
function evaluate_bicubicspline_interpolant end
function evaluate_two_dimensional_interpolant end
function plot_topography end
function plot_topography_with_interpolation_knots end
export evaluate_bicubicspline_interpolant, evaluate_two_dimensional_interpolant,
       plot_topography, plot_topography_with_interpolation_knots

# Note, empty routines for RBF interpolation are included and exported. They are extended
# in `ext/KernelInterpolationExt.jl` where their implementations are found.
function RBFInterpolation end
function RBFInterpolation1D end
function RBFInterpolation2D end
export RBFInterpolation, RBFInterpolation1D, RBFInterpolation2D

# Note, types and empty routines for Lavery spline interpolation are included and exported.
# They are extended in `ext/JuMPExt.jl` where their implementations are found.
# abstract type AbstractLaverySpline end
"""
    LaverySpline1D

One dimensional Lavery spline storage.

The attributes are:
- `x`: Vector of data points in x-direction (knot points)
- `y`: Vector of data values in y-direction
- `b`: Vector of slope coefficients at each knot point
- `weight`: Regularization weight for the slope coefficients (default: 1e-4)
- `integral_steps`: Number of integration steps for computation (default: 10)
"""
mutable struct LaverySpline1D{T <: Real}
    x::Vector{T}
    y::Vector{T}
    b::Vector{T}
    weight::T
    integral_steps::Int
end
"""
    LaverySpline2D

Two dimensional Lavery spline storage.

The attributes are:
- `x`: Vector of data points in x-direction (knot points)
- `y`: Vector of data points in x-direction (knot points)
- `z`: Matrix of data values in the z-direction
- `bx`: Matrix of slope coefficients at each knot point
- `by`: Matrix of slope coefficients at each knot point
- `lambda`: Additional regularization weight
"""
struct LaverySpline2D{T <: Real}
    x::Vector{T}
    y::Vector{T}
    z::Matrix{T}
    bx::Matrix{T}
    by::Matrix{T}
    lambda::T
end
function spline_interpolation(::Union{LaverySpline1D, LaverySpline2D}, args...)
    throw(MethodError(spline_interpolation, (args...,)))
end
function LaverySpline1D end
function LaverySpline2D end
export LaverySpline1D, LaverySpline2D
end
