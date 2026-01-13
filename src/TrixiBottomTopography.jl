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
using JuMP: Model, @variable, @objective, @constraint, optimize!, termination_status, OPTIMAL, value, set_silent, set_attribute
using HiGHS

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
export LinearBSpline, CubicBSpline, LaverySpline
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
end
