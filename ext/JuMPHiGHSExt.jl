# Package extension for adding JuMP modeling language and HiGHS optimization features
# for Lavery spline interpolation to TrixiBottomTopography.jl
module JuMPHiGHSExt

using JuMP: Model, VariableRef, direct_model, @variable, @objective, @constraint,
            optimize!, value, set_silent, set_attribute
using HiGHS

import MathOptInterface as MOI

using TrixiBottomTopography

import TrixiBottomTopography: LaverySpline1D, LaverySpline2D, spline_interpolation

##########################################
# One dimensional Lavery spline routines #
##########################################

# Helper container so that the JuMP model only needs constructed once
mutable struct LaverySpline1DModel{T <: Real}
    model::Model
    b::Vector{VariableRef}
    abs_b::Vector{VariableRef}
    abs_E::Matrix{VariableRef}
    delta_y::Vector{T} # data-dependent coefficients
end

# Constructor
function LaverySpline1DModel(len::Int, lambda::T, integral_steps::Int) where {T <: Real}
    sumDomain = 1:(len - 1)
    bDomain = 1:len

    # integration grid
    integralDomain = collect(range(convert(T, -0.5), convert(T, 0.5);
                                   length = integral_steps))

    # precompute t-dependent scalars
    # that are the tridiagonal matrix entries for a spline
    aa = convert(T, -1) .+ convert(T, 6) .* integralDomain
    bb = convert(T, 1) .+ convert(T, 6) .* integralDomain
    cc = convert(T, 12) .* integralDomain

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
    delta_y = zeros(T, len - 1)

    inv_steps = convert(T, 1) / integral_steps

    @objective(model, Min,
               sum(inv_steps * abs_E[i, k] for i in sumDomain, k in 1:integral_steps)+
               sum(lambda * abs_b[i] for i in bDomain))

    # constraints
    @constraint(model, [i in sumDomain, k in 1:integral_steps],
                abs_E[i, k]>=aa[k] * b[i] + bb[k] * b[i + 1] - cc[k] * delta_y[i])

    @constraint(model, [i in sumDomain, k in 1:integral_steps],
                abs_E[i, k]>=-aa[k] * b[i] - bb[k] * b[i + 1] + cc[k] * delta_y[i])

    @constraint(model, [i in bDomain], abs_b[i]>=b[i])
    @constraint(model, [i in bDomain], abs_b[i]>=-b[i])

    return LaverySpline1DModel{T}(model, b, abs_b, abs_E, delta_y)
end

"""
    LaverySpline1D(x::Vector, y::Vector; lambda::Real = 0.0, integral_steps::Int = 10)

This function calculates the inputs for the structure [`LaverySpline1D`](@ref).
The input values are:
- `x`: Vector of x-coordinates of the data points (knots)
- `y`: Vector of y-coordinates (function values) at the data points
- `lambda`: Additional regularization parameter for smoothness (default: 0.0)
- `integral_steps`: Number of discrete points for integration (default: 10)

First the data is sorted via [`sort_data`](@ref TrixiBottomTopography.sort_data)
to guarantee that the `x` and `y` values are in ascending order.

The Lavery spline coefficients are computed by solving an optimization problem.
The objective function in the optimization procedure is the second derivative of the spline.
to control the local curvature (point-wise).
Additional positivity constraints are placed on the coefficients to ensure no new extrema
are created by the resulting spline.

The optimization problem is formulated as a linear program and solved using HiGHS.jl.

References:
- John E. Lavery (2000),
   Univariate cubic Lp splines and shape-preserving, multiscale interpolation by univariate cubic L1 splines.
   [DOI: 10.1016/S0167-8396(00)00003-0](https://doi.org/10.1016/S0167-8396(00)00003-0)
- Logan Eriksson and Oscar Jemsson (2024)
   Uni- and bivariate interpolation of multiscale data using cubic L1 splines.
   [DiVA1918338](https://www.diva-portal.org/smash/get/diva2:1918338/FULLTEXT01.pdf)
"""
function LaverySpline1D(x::Vector{T}, y::Vector{T}; lambda::T = convert(T, 0.0),
                        integral_steps::Int = 10) where {T <: Real}
    if length(x) != length(y)
        throw(DimensionMismatch("Vectors x and y have to contain the same number of values"))
    end

    if length(x) < 2
        throw(ArgumentError("To perform Lavery spline interpolation, we need a vector
                             that contains at least 2 values."))
    end

    # Sort data to ensure ascending x values
    x, y = TrixiBottomTopography.sort_data(collect(x), collect(y))

    len = length(x)

    # Build the JuMP model once to save time
    spline_model = LaverySpline1DModel(len, lambda, integral_steps)

    for i in 1:(len - 1)
        hi = x[i + 1] - x[i]
        spline_model.delta_y[i] = (y[i + 1] - y[i]) / hi
    end

    # Solve optimization problem
    optimize!(spline_model.model)

    # Extract solution
    b_values = Vector{Float64}(undef, len)
    for i in 1:len
        b_values[i] = value(spline_model.b[i])
    end

    LaverySpline1D(x, y, b_values, lambda, integral_steps)
end

"""
    LaverySpline1D(path::String; lambda::Real = 0.0, integral_steps::Int = 10)

A function that reads in the `x` and `y` values for [`LaverySpline1D`](@ref)
from a .txt file.
The input values are:
- `path`: String of a path of the specific .txt file
- `lambda`: Additional regularization parameter for smoothness (default: 0.0)
- `integral_steps`: Number of discrete points for integration (default: 10)

The .txt file has to have the following structure to be interpreted by this function:
- First line: comment `# Number of x values`
- Second line: integer which gives the number of `x` values
- Third line: comment `# x values`
- Following lines: the `x` values where each value has its own line
- Line after the x-values: comment `# y values`
- Remaining lines: `y` values where each value has its own line

Note that the number of `x` and `y` values have to be the same.
An example can be found [here](https://gist.githubusercontent.com/maxbertrand1996/b05a90e66025ee1ebddf444a32c3fa01/raw/90d375c1ac11b26589aab1fe92bd0e6f6daf37b7/Rhine_data_1D_10.txt)
"""
function LaverySpline1D(path::String; kwargs...)
    x, y = TrixiBottomTopography.parse_txt_1D(path)

    LaverySpline1D(x, y; kwargs...)
end

"""
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
function spline_interpolation(lavery_spline::LaverySpline1D, x::Number)
    xData = lavery_spline.x
    yData = lavery_spline.y
    b = lavery_spline.b

    # Find the patch containing x
    i = max(1, min(searchsortedlast(xData, x), length(xData) - 1))

    # Helper functions for this patch
    h_i = xData[i + 1] - xData[i]
    dy_i = (yData[i + 1] - yData[i]) / h_i

    # Local coordinate within patch
    x_i = x - xData[i]

    # Evaluate cubic polynomial (written with Horner's rule)
    y_val = (yData[i] +
             x_i * (b[i] +
              x_i * ((-(2 * b[i] + b[i + 1]) + 3 * dy_i) / h_i +
               (b[i] + b[i + 1] - 2 * dy_i) * x_i / h_i^2)))

    return y_val
end

##########################################
# Two dimensional Lavery spline routines #
##########################################

"""
    LaverySpline2D(x::Vector, y::Vector, z::Matrix; lambda::Real = 0.0)

This function calculates the inputs for the structure [`LaverySpline2D`](@ref LaverySpline2D).
The input values are:
- `x`: Vector of x-coordinates of the data points (knots)
- `y`: Vector of y-coordinates of the data points (knots)
- `z`: Matrix of z-coordinates (function values) at the data points
- `lambda`: Additional regularization parameter for smoothness (default: 0.0)

Creates a total variation (TV) regularized bicubic spline.
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

The data is sorted via [`sort_data`](@ref TrixiBottomTopography.sort_data)
to guarantee that the `x`, `y`, and `z` values are in ascending order.

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
function LaverySpline2D(x::Vector{T}, y::Vector{T}, z::Matrix{T};
                        lambda::T = convert(T, 0)) where {T <: Real}
    n = length(x)
    m = length(y)

    if size(z, 2) != n || size(z, 1) != m
        error("Dimension mismatch between x, y, z")
    end

    x, y, z = TrixiBottomTopography.sort_data(x, y, z)

    # As with the `BicubicBSpline` the data is passed in the ordering
    # z[j, i] where j = 1, ..., m and i = 1, ..., n
    # So we transpose into the format needed by the implementation below
    z = Matrix{T}(transpose(z))

    # offsets
    bx_off = 0
    by_off = n * m
    dxbx_off = 2 * n * m
    dyby_off = dxbx_off + (n-1) * m
    dybx_off = dyby_off + n * (m-1)
    dxby_off = dybx_off + n * (m-1)
    n_vars = dxby_off + (n-1) * m
        
    @inline bx_idx(i,j) = bx_off + (i-1)*m + j
    @inline by_idx(i,j) = by_off + (i-1)*m + j
    @inline dxbx_idx(i, j) = dxbx_off + (i - 1) * m + j
    @inline dyby_idx(i, j) = dyby_off + (i - 1) * (m - 1) + j 
    @inline dybx_idx(i,j) = dybx_off + (i-1) * (m-1) + j
    @inline dxby_idx(i,j) = dxby_off + (i-1) * m + j

    # Setup the model and solver choice
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, MOI.RawOptimizerAttribute("solver"),  "simplex")
    MOI.set(model, MOI.RawOptimizerAttribute("presolve"), "off")

    # Add variables
    vars = MOI.add_variables(model, n_vars)

    lb_zero = MOI.GreaterThan(zero(T))
    for v in vars
        MOI.add_constraint(model, v, lb_zero)
    end

    # Coefficients 
    obj_terms = Vector{MOI.ScalarAffineTerm{T}}(undef, n_vars)

    for k in 1:n_vars
        c = if k <= 2 * n * m
            lambda            # bx or by block
        else
            one(T)            # dxbx, dyby, dyxbx, dxyby blocks
        end
        obj_terms[k] = MOI.ScalarAffineTerm{T}(c, vars[k])
    end

    # First-order total variation (TV) of gradients
    # Additional regularization can be added with `lambda > 0`
    MOI.set(model,
    MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
    MOI.ScalarAffineFunction{T}(obj_terms, zero(T)))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    
    # Constraints
    nnz_per_con = 3
    buf = Vector{MOI.ScalarAffineTerm{T}}(undef, nnz_per_con)
    rhs = MOI.GreaterThan(zero(T))

    # abs_dxbx: s_{i,j} ≥ ± (bx[i+1,j] - bx[i,j]) for i = 1:n-1, j = 1:m
    for i in 1:(n-1), j in 1:m
        s = vars[dxbx_idx(i, j)]
        bp = vars[bx_idx(i+1, j)]
        bm = vars[bx_idx(i, j)]
        # s - bp + bm ≥ 0
        buf[1] = MOI.ScalarAffineTerm{T}(one(T),s)
        buf[2] = MOI.ScalarAffineTerm{T}(-one(T),bp)
        buf[3] = MOI.ScalarAffineTerm{T}(one(T),bm)
        MOI.add_constraint(model, MOI.ScalarAffineFunction{T}(copy(buf), zero(T)), rhs)
        # s + bp - bm ≥ 0
        buf[1] = MOI.ScalarAffineTerm{T}(one(T),s)
        buf[2] = MOI.ScalarAffineTerm{T}(one(T),bp)
        buf[3] = MOI.ScalarAffineTerm{T}(-one(T),bm)
        MOI.add_constraint(model, MOI.ScalarAffineFunction{T}(copy(buf), zero(T)), rhs)
    end 

    # abs_dyby: s_{i,j} ≥ ± (by[i,j+1] - by[i,j]) for i = 1:n, j = 1:m-1
    for i in 1:n, j in 1:(m-1)
        s = vars[dyby_idx(i,j)]
        bp = vars[by_idx(i, j+1)]
        bm = vars[by_idx(i,j)]
        # s - bp + bm ≥ 0
        buf[1] = MOI.ScalarAffineTerm{T}(one(T),s)
        buf[2] = MOI.ScalarAffineTerm{T}(-one(T),bp)
        buf[3] = MOI.ScalarAffineTerm{T}(one(T),bm)
        MOI.add_constraint(model, MOI.ScalarAffineFunction{T}(copy(buf), zero(T)), rhs)
        # s + bp - bm ≥ 0
        buf[1] = MOI.ScalarAffineTerm{T}(one(T),s)
        buf[2] = MOI.ScalarAffineTerm{T}(one(T),bp)
        buf[3] = MOI.ScalarAffineTerm{T}(-one(T),bm)
        MOI.add_constraint(model, MOI.ScalarAffineFunction{T}(copy(buf), zero(T)), rhs)
    end

    # abs_dybx: s_{i,j} ≥ ± (bx[i,j+1] - bx[i,j]) for i = 1:n, j = 1:m-1
    for i in 1:n, j in 1:(m-1)
        s = vars[dybx_idx(i,j)]
        bp = vars[bx_idx(i,j+1)]
        bm = vars[bx_idx(i,j)]
        # s - bp + bm ≥ 0
        buf[1] = MOI.ScalarAffineTerm{T}(one(T),s)
        buf[2] = MOI.ScalarAffineTerm{T}(-one(T),bp)
        buf[3] = MOI.ScalarAffineTerm{T}(one(T), bm)
        MOI.add_constraint(model, MOI.ScalarAffineFunction{T}(copy(buf), zero(T)), rhs)
        # s + bp - bm ≥ 0
        buf[1] = MOI.ScalarAffineTerm{T}(one(T),s)
        buf[2] = MOI.ScalarAffineTerm{T}(one(T),bp)
        buf[3] = MOI.ScalarAffineTerm{T}(-one(T), bm)
        MOI.add_constraint(model, MOI.ScalarAffineFunction{T}(copy(buf), zero(T)), rhs)
    end

    # abs_dxby: s_{i,j} ≥ ± (by[i+1,j] - by[i,j]) for i = 1:n-1, j = 1:m
    for i in 1:(n-1), j in 1:m
        s = vars[dxby_idx(i,j)]
        bp = vars[by_idx(i+1, j)]
        bm = vars[by_idx(i, j)]
        # s - bp + bm ≥ 0
        buf[1] = MOI.ScalarAffineTerm{T}(one(T),s)
        buf[2] = MOI.ScalarAffineTerm{T}(-one(T), bp)
        buf[3] = MOI.ScalarAffineTerm{T}(one(T), bm)
        MOI.add_constraint(model, MOI.ScalarAffineFunction{T}(copy(buf), zero(T)), rhs)
        # s + bp - bm ≥ 0
        buf[1] = MOI.ScalarAffineTerm{T}(one(T),s)
        buf[2] = MOI.ScalarAffineTerm{T}(one(T), bp)
        buf[3] = MOI.ScalarAffineTerm{T}(-one(T), bm)
        MOI.add_constraint(model, MOI.ScalarAffineFunction{T}(copy(buf), zero(T)), rhs)
    end

    # Solve
    MOI.optimize!(model)

    # Extract solution
    bxv = Matrix{T}(undef, n, m)
    byv = Matrix{T}(undef, n, m)
    for i in 1:n, j in 1:m
        bxv[i, j] = MOI.get(model, MOI.VariablePrimal(), vars[bx_idx(i, j)])
        byv[i, j] = MOI.get(model, MOI.VariablePrimal(), vars[by_idx(i, j)])
    end

    return LaverySpline2D(x, y, z, bxv, byv, lambda)
end

"""
    LaverySpline2D(path::String; lambda::Real = 0.0)

A function which reads in the `x`, `y` and `z` values for [`LaverySpline2D`](@ref)
from a .txt file.
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
function LaverySpline2D(path::String; kwargs...)
    x, y, z = TrixiBottomTopography.parse_txt_2D(path)

    LaverySpline2D(x, y, z; kwargs...)
end

@doc raw"""
    spline_interpolation(lavery_spline::LaverySpline2D, x::Number, y::Number)

Evaluates a total variation (TV) based spline at a single point `(x, y)`.
The TV regularized spline behaves as a bicubic Lavery spline in that no new extrema are generated
and shape is preserved.
The coefficients are computed through an optimization procedure,
see [`LaverySpline2D`](@ref).
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
function spline_interpolation(lavery_spline::LaverySpline2D, x::Number, y::Number)
    xData = lavery_spline.x
    yData = lavery_spline.y
    z = lavery_spline.z
    bx = lavery_spline.bx
    by = lavery_spline.by

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
