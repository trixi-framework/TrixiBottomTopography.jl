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
include("1D/spline_cache_1D.jl")
include("1D/spline_methods_1D.jl")
include("1D/spline_utils_1D.jl")

# Include two dimensional B-spline interpolation
include("2D/spline_cache_2D.jl")
include("2D/spline_methods_2D.jl")
include("2D/spline_utils_2D.jl")

# Include auxiliary functions
include("auxiliary/convert.jl")
include("auxiliary/default_example.jl")

# Export the functions which are used for B-spline interpolation
export LinearBSpline, CubicBSpline
export BilinearBSpline, BicubicBSpline
export spline_interpolation

# Export the functions which are used DGM and GeophysicalModelgenerator data conversion
export convert_dgm_1d, convert_dgm_2d, convert_geo_1d, convert_geo_2d

# Note, Empty routines for GeophysicalModelGenerator are included and exported. 
# They are extended in `ext/TrixiBottomTopographyGeophysicalModelGeneratorExt.jl`
function geo_topo_impression end
function create_topography_data end

# Note, Empty routines for visualization are included and exported. They are extended
# in `ext/TrixiBottomTopographyMakieExt.jl` where their implementations are found.
function evaluate_bicubicspline_interpolant end
function plot_topography end
function plot_topography_with_interpolation_knots end
export geo_topo_impression, create_topography_data
export evaluate_bicubicspline_interpolant,
       plot_topography, plot_topography_with_interpolation_knots

end