# Package extension for adding KernelInterpolation.jl-based features to TrixiBottomTopography.jl
module KernelInterpolationExt

using KernelInterpolation: NodeSet, interpolate
using TrixiBottomTopography

"""
    RBFInterpolation(nodeset::NodeSet, z, args...; kwargs...)

Create a radial basis function (RBF) interpolant using [KernelInterpolation.jl](https://github.com/JoshuaLampert/KernelInterpolation.jl)
from a set of scattered nodes and corresponding elevation values `z`.
The arguments `args` and keyword arguments `kwargs` are passed to the `interpolate` function of KernelInterpolation.jl.
See the [KernelInterpolation.jl documentation](https://joshualampert.github.io/KernelInterpolation.jl/stable/) for more details.
"""
function TrixiBottomTopography.RBFInterpolation(nodeset::NodeSet, z, args...; kwargs...)
    return interpolate(nodeset, z, args...; kwargs...)
end

"""
    RBFInterpolation1D(x, y, args...; kwargs...)
    RBFInterpolation1D(path, args...; kwargs...)

Create a 1D radial basis function (RBF) interpolant using [KernelInterpolation.jl](https://github.com/JoshuaLampert/KernelInterpolation.jl).
The arguments `args` and keyword arguments `kwargs` are passed to the `interpolate` function of KernelInterpolation.jl.
See the [KernelInterpolation.jl documentation](https://joshualampert.github.io/KernelInterpolation.jl/stable/) for more details.

The RBF interpolant can be created with [`RBFInterpolation`](@ref) from a set of scattered nodes and corresponding elevation values `z` using
a `NodeSet` from KernelInterpolation.jl. Alternatively, the same interface as for [`BicubicBSpline`](@ref) and
[`BilinearBSpline`](@ref) can be used with `RBFInterpolation1D`, i.e., either by providing vectors `x` and `y` of coordinates and a matrix `z`
of elevation values or by providing a path to a text file containing the data `x`, `y`, and `z`.
"""
function TrixiBottomTopography.RBFInterpolation1D(x, y, args...; kwargs...)
    nodeset = NodeSet(vec([[xx] for xx in x]))
    return TrixiBottomTopography.RBFInterpolation(nodeset, vec(y), args...; kwargs...)
end

function TrixiBottomTopography.RBFInterpolation1D(path::String, args...; kwargs...)
    x, y = TrixiBottomTopography.parse_txt_1D(path)
    return TrixiBottomTopography.RBFInterpolation1D(x, y, args...; kwargs...)
end

"""
    RBFInterpolation2D(x, y, z, args...; kwargs...)
    RBFInterpolation2D(path, args...; kwargs...)

Create a 2D radial basis function (RBF) interpolant using [KernelInterpolation.jl](https://github.com/JoshuaLampert/KernelInterpolation.jl).
The arguments `args` and keyword arguments `kwargs` are passed to the `interpolate` function of KernelInterpolation.jl.
See the [KernelInterpolation.jl documentation](https://joshualampert.github.io/KernelInterpolation.jl/stable/) for more details.

The RBF interpolant can be created with [`RBFInterpolation`](@ref) from a set of scattered nodes and corresponding elevation values `z` using
a `NodeSet` from KernelInterpolation.jl. Alternatively, the same interface as for [`BicubicBSpline`](@ref) and
[`BilinearBSpline`](@ref) can be used with `RBFInterpolation2D`, i.e., either by providing vectors `x` and `y` of coordinates and a matrix `z`
of elevation values or by providing a path to a text file containing the data `x`, `y`, and `z`.
"""
function TrixiBottomTopography.RBFInterpolation2D(x, y, z, args...; kwargs...)
    nodeset = NodeSet(vec([[xx; yy] for xx in x, yy in y]))
    return TrixiBottomTopography.RBFInterpolation(nodeset, vec(Matrix(z')), args...;
                                                  kwargs...)
end

function TrixiBottomTopography.RBFInterpolation2D(path::String, args...; kwargs...)
    x, y, z = TrixiBottomTopography.parse_txt_2D(path)
    return TrixiBottomTopography.RBFInterpolation2D(x, y, z, args...; kwargs...)
end

end
