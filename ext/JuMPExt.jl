# Package extension for adding JuMP modeling language and HiGHS optimization features
# for Lavery spline interpolation to TrixiBottomTopography.jl
module JuMPExt

using JuMP: Model, VariableRef, direct_model, @variable, @objective, @constraint,
            optimize!, value, set_silent, set_attribute
using HiGHS

# Use exported function names `LaverySpline1D`, `LaverySpline2D`, and `spline_interpolation`
# to avoid namespace conflicts
using TrixiBottomTopography

import TrixiBottomTopography: AbstractLaverySpline

##########################################
# One dimensional Lavery spline routines #
##########################################

# Helper container so that the JuMP model only needs constructed once
mutable struct LaverySpline1DModel
    model::Model
    b::Vector{VariableRef}
    abs_b::Vector{VariableRef}
    abs_E::Matrix{VariableRef}
    deltaZ::Vector{Float64} # data-dependent coefficients
end

# Constructor
function build_lavery_spline_1d_model(len::Int, weight::Float64, integral_steps::Int)
    sumDomain = 1:(len - 1)
    bDomain = 1:len

    # integration grid
    integralDomain = collect(range(-0.5, 0.5; length = integral_steps))

    # precompute t-dependent scalars
    aa = -1 .+ 6 .* integralDomain
    bb = 1 .+ 6 .* integralDomain
    cc = 12 .* integralDomain

    # Setup the model and solver choice
    model = direct_model(HiGHS.Optimizer())
    set_attribute(model, "solver", "simplex")
    set_attribute(model, "presolve", "off")
    set_silent(model)

    # variables
    @variable(model, b[bDomain])
    @variable(model, abs_b[bDomain]>=0)
    @variable(model, abs_E[i in sumDomain, k in 1:integral_steps]>=0)

    # placeholder for data
    deltaZ = zeros(len - 1)

    inv_steps = 1 / integral_steps

    @objective(model, Min,
               sum(inv_steps * abs_E[i, k] for i in sumDomain, k in 1:integral_steps)+
               sum(weight * abs_b[i] for i in bDomain))

    # constraints
    @constraint(model, [i in sumDomain, k in 1:integral_steps],
                abs_E[i, k]>=aa[k] * b[i] + bb[k] * b[i + 1] - cc[k] * deltaZ[i])

    @constraint(model, [i in sumDomain, k in 1:integral_steps],
                abs_E[i, k]>=-aa[k] * b[i] - bb[k] * b[i + 1] + cc[k] * deltaZ[i])

    @constraint(model, [i in bDomain], abs_b[i]>=b[i])
    @constraint(model, [i in bDomain], abs_b[i]>=-b[i])

    return LaverySpline1DModel(model,
                               b, abs_b, abs_E,
                               deltaZ)
end

"""
    LaverySpline1D(x, z, b, weight, integral_steps)

One dimensional Lavery spline structure.

The attributes are:
- `x`: Vector of data points in x-direction (knot points)
- `z`: Vector of data values in z-direction
- `b`: Vector of slope coefficients at each knot point
- `weight`: Regularization weight for the slope coefficients (default: 1e-4)
- `integral_steps`: Number of integration steps for computation (default: 10)

"""
mutable struct LaverySpline1D{x_type, z_type, b_type} <: AbstractLaverySpline
    x::x_type
    z::z_type
    b::b_type
    weight::Float64
    integral_steps::Int
end

@doc raw"""
    LaverySpline1D(xData::AbstractVector, zData::AbstractVector;
                   weight::Float64 = 1e-4,
                   integral_steps::Int = 10)

This function calculates the inputs for the structure [`LaverySpline1D`](@ref).
The input values are:
- `xData`: Vector of x-coordinates of the data points (knots)
- `zData`: Vector of z-coordinates (function values) at the data points
- `weight`: Regularization parameter for smoothness (default: 1e-4)
- `integral_steps`: Number of discrete points for integration (default: 10)

First the data is sorted via [`sort_data`](@ref) to guarantee that the `x` values
are in ascending order.

The Lavery spline is computed by solving an optimization problem.

The optimization problem is formulated as a linear program and solved using HiGHS.jl.

References:
- John E. Lavery (2000),
   Univariate cubic Lp splines and shape-preserving, multiscale interpolation by univariate cubic L1 splines.
   [DOI: 10.1016/S0167-8396(00)00003-0](https://doi.org/10.1016/S0167-8396(00)00003-0)
- Logan Eriksson and Oscar Jemsson (2024)
   Uni- and bivariate interpolation of multiscale data using cubic L1 splines.
   [DiVA1918338](https://www.diva-portal.org/smash/get/diva2:1918338/FULLTEXT01.pdf)
"""
function TrixiBottomTopography.LaverySpline1D(xData::AbstractVector, zData::AbstractVector;
                                              weight::Float64 = 1e-4,
                                              integral_steps::Int = 10)
    if length(xData) != length(zData)
        throw(DimensionMismatch("Vectors xData and zData have to contain the same number of values"))
    end

    if length(xData) < 2
        throw(ArgumentError("To perform Lavery spline interpolation, we need an vector
                             that contains at least 2 values."))
    end

    # Sort data to ensure ascending x values
    xData, zData = TrixiBottomTopography.sort_data(collect(xData), collect(zData))

    len = length(xData)

    # Build the JuMP model once to save time
    spline_model = build_lavery_spline_1d_model(len, weight, integral_steps)

    for i in 1:(len - 1)
        hi = xData[i + 1] - xData[i]
        spline_model.deltaZ[i] = (zData[i + 1] - zData[i]) / hi
    end

    # Solve optimization problem
    optimize!(spline_model.model)

    # Extract solution
    b_values = Vector{Float64}(undef, len)
    for i in 1:len
        b_values[i] = value(spline_model.b[i])
    end

    LaverySpline1D(xData, zData, b_values, weight, integral_steps)
end

"""
    LaverySpline1D(path::String; weight::Float64 = 1e-4, integral_steps::Int = 10)

A function that reads in the `x` and `z` values for [`LaverySpline1D`](@ref) from a .txt file.
The input values are:
- `path`: String of a path of the specific .txt file
- `weight`: Regularization parameter for smoothness (default: 1e-4)
- `integral_steps`: Number of discrete points for integration (default: 10)

The .txt file has to have the following structure to be interpreted by this function:
- First line: comment `# Number of x values`
- Second line: integer which gives the number of `x` values
- Third line: comment `# x values`
- Following lines: the `x` values where each value has its own line
- Line after the x-values: comment `# z values`
- Remaining lines: `z` values where each value has its own line

Note that the number of `x` and `z` values have to be the same.
An example can be found [here](https://gist.githubusercontent.com/maxbertrand1996/b05a90e66025ee1ebddf444a32c3fa01/raw/90d375c1ac11b26589aab1fe92bd0e6f6daf37b7/Rhine_data_1D_10.txt)
"""
function TrixiBottomTopography.LaverySpline1D(path::String; kwargs...)
    x, z = TrixiBottomTopography.parse_txt_1D(path)

    TrixiBottomTopography.LaverySpline1D(x, z; kwargs...)
end

@doc raw"""
    spline_interpolation(lavery_spline::LaverySpline1D, x::Number)

Evaluate the one dimensional Lavery spline at a single point `x`.
The form of the spline originally stated in Lavery is given in equation (2.3a)
of the thesis of Eriksson and Jemsson.

References:
- John E. Lavery (2000),
   Univariate cubic Lp splines and shape-preserving, multiscale interpolation by univariate cubic L1 splines.
   [DOI: 10.1016/S0167-8396(00)00003-0](https://doi.org/10.1016/S0167-8396(00)00003-0)
- Logan Eriksson and Oscar Jemsson (2024)
   Uni- and bivariate interpolation of multiscale data using cubic L1 splines.
   [DiVA1918338](https://www.diva-portal.org/smash/get/diva2:1918338/FULLTEXT01.pdf)
"""
function TrixiBottomTopography.spline_interpolation(lavery_spline::LaverySpline1D,
                                                    x::Number)
    xData = lavery_spline.x
    zData = lavery_spline.z
    b = lavery_spline.b

    # Find the patch containing x
    i = max(1, min(searchsortedlast(xData, x), length(xData) - 1))

    # Helper functions for this patch
    h_i = xData[i + 1] - xData[i]
    dz_i = (zData[i + 1] - zData[i]) / h_i

    # Local coordinate within patch
    x_i = x - xData[i]

    # Evaluate cubic polynomial (written with Horner's rule)
    z_val = (zData[i] +
             x_i * (b[i] +
              x_i * ((-(2 * b[i] + b[i + 1]) + 3 * dz_i) / h_i +
               (b[i] + b[i + 1] - 2 * dz_i) * x_i / h_i^2)))

    return z_val
end

##########################################
# Two dimensional Lavery spline routines #
##########################################

"""
    LaverySpline2D(x, y, z, bx, by, lambda)

Two dimensional Lavery spline structure.

The attributes are:
- `x`: Vector of data points in x-direction (knot points)
- `y`: Vector of data points in x-direction (knot points)
- `z`: Matrix of data values in the z-direction
- `bx`: Matrix of slope coefficients at each knot point
- `by`: Matrix of slope coefficients at each knot point
- `lambda`: Additional regularization weight
"""
struct LaverySpline2D <: AbstractLaverySpline
    x::AbstractVector{Float64}
    y::AbstractVector{Float64}
    z::AbstractMatrix{Float64}
    bx::AbstractMatrix{Float64}
    by::AbstractMatrix{Float64}
    lambda::Float64
end

@doc raw"""
    LaverySpline2D(xData::AbstractVector, yData::AbstractVector,
                   zData::AbstractMatrix;
                   lambda::Float64 = 0.0)

This function calculates the inputs for the structure [`LaverySpline2D`](@ref).
The input values are:
- `xData`: Vector of x-coordinates of the data points (knots)
- `yData`: Vector of y-coordinates of the data points (knots)
- `zData`: Matrix of z-coordinates (function values) at the data points
- `lambda`: Additional regularization parameter for smoothness (default: 0.0)

Creates a total variation (TV) bicubic spline.
The spline is written as a tensor product of cubic Hermite basis functions.
The coefficients for the spline are computed from an optimization problem with
constraints to ensure no new extrema and shape preservation.

The original formulation from Lavery uses an objective function in terms of the second derivatives
to control curvature of the resulting spline pointwise.
This strategy requires geometric subdivision via Sibson elements and is numerically expensive.
This implementation, instead, uses sign constraints on the gradients together with a TV penalty
and an optimization procedure to determine the coefficients and maintain shape preservation globally.
In essence, this is a discrete convex relaxation of the original
formulation of Lavery that is tight under monotonicity.

The data is sorted via [`sort_data`](@ref) to guarantee that the `x`, `y`, and `z`
values are in ascending order.

The optimization problem is formulated as a linear program and solved using
the simplex method available in HiGHS.jl.

References:
- John E. Lavery (2001),
   Shape-preserving, multiscale interpolation by bi- and multivariate cubic splines.
   [DOI: 10.1016/S0167-8396(01)00034-6](https://doi.org/10.1016/S0167-8396(01)00034-6)
- Logan Eriksson and Oscar Jemsson (2024)
   Uni- and bivariate interpolation of multiscale data using cubic L1 splines.
   [DiVA1918338](https://www.diva-portal.org/smash/get/diva2:1918338/FULLTEXT01.pdf)
- Frans Kuijt and Ruud van Damme (2001)
   A linear approach to shape preserving spline approximation
   [DOI: 10.1023/A:1016660513725](https://doi.org/10.1023/A:1016660513725)
- Lu Yu, Qingwei Jin, John E. Lavery, and Shu-Cherng Fang (2010)
   Univariate Cubic L1 Interpolating Splines: Spline Functional, Window Size and Analysis-based Algorithm.
   [DOI: 10.3390/a3030311](https://doi.org/10.3390/a3030311)
- Ziteng Wang, John Lavery, and Shu-Cherng Fang (2014)
   Approximation of Irregular Geometric Data by Locally Calculated Univariate Cubic L1 Spline Fits
   [DOI: 10.1007/s40745-014-0002-z](https://doi.org/10.1007/s40745-014-0002-z)
"""
function TrixiBottomTopography.LaverySpline2D(xData::AbstractVector, yData::AbstractVector,
                                              zData::AbstractMatrix; lambda::Float64 = 0.0)
    # TODO: maybe change to nx and ny or something or n and m
    I = length(xData)
    J = length(yData)

    if size(zData, 2) != I || size(zData, 1) != J
        error("Dimension mismatch between x, y, z")
    end

    xData, yData, zData = TrixiBottomTopography.sort_data(xData, yData, zData)

    # As with the `BicubicBSpline` the data is passed in the ordering
    # zData[j, i] where j = 1, ..., J and i = 1, ..., I
    # So we transpose into the format needed by the implementation below
    zData = transpose(zData)

    # Setup the model and solver choice
    model = direct_model(HiGHS.Optimizer())
    set_attribute(model, "solver", "simplex")
    set_attribute(model, "presolve", "off")
    set_silent(model)

    # Setup the variables that enforce monotonicity via bounds
    @variable(model, bx[1:I, 1:J]>=0)
    @variable(model, by[1:I, 1:J]>=0)

    @variable(model, abs_dxbx[1:(I - 1), 1:J]>=0)
    @variable(model, abs_dyby[1:I, 1:(J - 1)]>=0)
    @variable(model, abs_dybx[1:I, 1:(J - 1)]>=0)
    @variable(model, abs_dxby[1:(I - 1), 1:J]>=0)

    # Absolute value constraints on the gradients
    @constraint(model, [i = 1:(I - 1), j = 1:J], abs_dxbx[i, j]>=bx[i + 1, j] - bx[i, j])
    @constraint(model, [i = 1:(I - 1), j = 1:J], abs_dxbx[i, j]>=-bx[i + 1, j] + bx[i, j])

    @constraint(model, [i = 1:I, j = 1:(J - 1)], abs_dyby[i, j]>=by[i, j + 1] - by[i, j])
    @constraint(model, [i = 1:I, j = 1:(J - 1)], abs_dyby[i, j]>=-by[i, j + 1] + by[i, j])

    @constraint(model, [i = 1:I, j = 1:(J - 1)], abs_dybx[i, j]>=bx[i, j + 1] - bx[i, j])
    @constraint(model, [i = 1:I, j = 1:(J - 1)], abs_dybx[i, j]>=-bx[i, j + 1] + bx[i, j])

    @constraint(model, [i = 1:(I - 1), j = 1:J], abs_dxby[i, j]>=by[i + 1, j] - by[i, j])
    @constraint(model, [i = 1:(I - 1), j = 1:J], abs_dxby[i, j]>=-by[i + 1, j] + by[i, j])

    # Objective function is the first-order total variation (TV) of gradients
    # Additional regularization can be added with `lambda > 0`
    @objective(model, Min,
               sum(abs_dxbx)+sum(abs_dyby)
               +sum(abs_dybx)+sum(abs_dxby)
               +lambda*(sum(bx) + sum(by)))

    # Solve optimization problem
    optimize!(model)

    # Extract solution
    bxv = value.(bx)
    byv = value.(by)

    return LaverySpline2D(xData, yData, zData, bxv, byv, lambda)
end

"""
    LaverySpline2D(path::String; lambda::Float64 = 0.0)

A function which reads in the `x`, `y` and `z` values for [`LaverySpline2D`](@ref) from a .txt file.
The input values are:
- `path`: String of a path of the specific .txt file
- `lambda`: Regularization parameter for smoothness (default: 0.0)

The .txt file has to have the following structure to be interpreted by this function:
- First line: comment `# Number of x values`
- Second line: integer which gives the number of `x` values
- Third line: comment `# Number of y values`
- Fourth line: integer which gives the number of `y` values
- Fifth line: comment `# x values`
- Following lines: the `x` values where each value has its own line
- Line after the x-values: comment `# y values`
- Following lines: `y` values where each value has its own line
- Line after the y-values: comment `# z values`
- Remaining lines: values for `z` where each value has its own line and is in th following order:
                   z_11, z_12, ... z_1n, z_21, ... z_2n, ..., z_m1, ..., z_mn
"""
function TrixiBottomTopography.LaverySpline2D(path::String; kwargs...)
    x, y, z = TrixiBottomTopography.parse_txt_2D(path)

    TrixiBottomTopography.LaverySpline2D(x, y, z; kwargs...)
end

@doc raw"""
    spline_interpolation(lavery_spline::LaverySpline2D, x::Number, y::Number)

Evaluates a total variation (TV) based spline at a single point `(x,y)`.
The TV spline behaves as a bicubic Lavery spline in that no new extrema are generated
and shape is preserved.
The coefficients are computed through an optimization procedure, see [`LaverySpline2D`](@ref).
The spline is formulated as a tensor product of cubic Hermite splines with basis functions
in each direction
```math
\begin{aligned}
   h_{00}(t) &= 2 t^3 - 3 t^2 + 1\\
   h_{10}(t) &= t^3 - 2 t^2 + t\\
   h_{01}(t) &= -2 t^3 + 3 t^2\\
   h_{11}(t) &= t^3 - t^2\\
\end{aligned}
```

References:
- Lu Yu, Qingwei Jin, John E. Lavery, and Shu-Cherng Fang (2010)
   Univariate Cubic L1 Interpolating Splines: Spline Functional, Window Size and Analysis-based Algorithm.
   [DOI: 10.3390/a3030311](https://doi.org/10.3390/a3030311)
- Ziteng Wang, John Lavery & Shu-Cherng Fang (2014)
   Approximation of Irregular Geometric Data by Locally Calculated Univariate Cubic L1 Spline Fits
   [DOI: 10.1007/s40745-014-0002-z](https://doi.org/10.1007/s40745-014-0002-z)
"""
function TrixiBottomTopography.spline_interpolation(lavery::LaverySpline2D, x::Number,
                                                    y::Number)
    xData = lavery.x
    yData = lavery.y
    z = lavery.z
    bx = lavery.bx
    by = lavery.by

    # Find containing cell
    i = max(1, min(searchsortedlast(xData, x), length(xData) - 1))
    j = max(1, min(searchsortedlast(yData, y), length(yData) - 1))

    dx = xData[i + 1] - xData[i]
    dy = yData[j + 1] - yData[j]

    xt = (x - xData[i]) / dx
    yt = (y - yData[j]) / dy

    # cubic Hermite spline basis functions
    h00(t) = 1 + t^2 * (-3 + 2 * t)
    h10(t) = t + t^2 * (-2 + t)
    h01(t) = t^2 * (3 - 2 * t)
    h11(t) = t^2 * (-1 + t)

    # Corner values of the current cell
    z00 = z[i, j]
    z10 = z[i + 1, j]
    z01 = z[i, j + 1]
    z11 = z[i + 1, j + 1]

    bx00 = bx[i, j]
    bx10 = bx[i + 1, j]
    bx01 = bx[i, j + 1]
    bx11 = bx[i + 1, j + 1]

    by00 = by[i, j]
    by10 = by[i + 1, j]
    by01 = by[i, j + 1]
    by11 = by[i + 1, j + 1]

    # Tensor product cubic Hermite interpolation
    return (h00(xt) * (h00(yt) * z00 + h01(yt) * z01) +
            h01(xt) * (h00(yt) * z10 + h01(yt) * z11) +
            dx * (h10(xt) * (h00(yt) * bx00 + h01(yt) * bx01) +
             h11(xt) * (h00(yt) * bx10 + h01(yt) * bx11)) +
            dy * (h00(xt) * (h10(yt) * by00 + h11(yt) * by01) +
             h01(xt) * (h10(yt) * by10 + h11(yt) * by11)))
end

end
