
#######################
### Linear B Spline ###
#######################

# Linear B Spline structure
"""
    LinearBSpline(x, Delta, Q, IP)

One dimensional linear B-spline structure which contains all important attributes to define
a B-Spline interpolation function.
These attributes are:
- `x`: Vector of values in x-direction
- `Delta`: Length of a single patch in the given data set. A patch is the area between two
           consecutive `x` values. The value `Delta` corresponds to the distance between two
           consecutive values in x-direction. As we are only considering Cartesian grids,
           `Delta` is equal for all patches
- `Q`: Vector which contains the control points
- `IP`: Coefficients matrix
"""
mutable struct LinearBSpline{x_type, Delta_type, Q_type, IP_type}
    x::x_type
    Delta::Delta_type
    Q::Q_type
    IP::IP_type
end

# Fill structure
@doc raw"""
    LinearBSpline(x::Vector, y::Vector)

This function calculates the inputs for the structure [`LinearBSpline`](@ref).
The input values are:
- `x`: A vector that contains equally spaced values in x-direction
- `y`: A vector that contains values in y-direction

Linear B-spline interpolation is only possible if the data set has at least two values in `x`.

First the data is sorted via [`sort_data`](@ref) to
guarantee that the `x` values are in ascending order.

The patch size `Delta` is calculated by subtracting the second and first `x` values.
This can be done because we only consider equally spaced `x` values.
A patch is the area between two consecutive `x` values.

For linear B-spline interpolation, the control points `Q` correspond with the values in `y`.

The coefficients matrix `IP` for linear B-splines is fixed to be
```math
\begin{aligned}
  \begin{pmatrix}
    -1 & 1\\
    1 & 0
  \end{pmatrix}
\end{aligned}
```

A reference for the calculations in this script can be found in Chapter 1 of
-  Quentin Agrapart & Alain Batailly (2020),
   Cubic and bicubic spline interpolation in Python.
   [hal-03017566v2](https://hal.archives-ouvertes.fr/hal-03017566v2)
"""
function LinearBSpline(x::Vector, y::Vector)
    if length(x) != length(y)
        throw(DimensionMismatch("Vectors x and y have to contain the same number of values"))
    end

    if length(x) == 1
        throw(ArgumentError("To perform linear B-spline interpolation, we need an x vector
                             that contains at least 2 values."))
    end

    x, y = sort_data(x, y)

    Delta = x[2] - x[1]

    IP = @SMatrix [-1 1;
                   1 0]

    Q = y

    LinearBSpline(x, Delta, Q, IP)
end

# Read from file
"""
    LinearBSpline(path::String)

A function that reads in the `x` and `y` values for
[`LinearBSpline`](@ref) from a .txt file.
The input values are:
- `path`: String of a path of the specific .txt file

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
function LinearBSpline(path::String)
    x, y = parse_txt_1D(path)

    LinearBSpline(x, y)
end

######################
### Cubic B Spline ###
######################

# Cubic B Spline structure
"""
    CubicBSpline(x, Delta, Q, IP)

One dimensional cubic B-spline structure that contains all important attributes to define
a B-Spline interpolation function. Similar to [`LinearBSpline`](@ref)
These attributes are:
- `x`: Vector of values in x-direction
- `Delta`: Length of a single patch in the given data set. A patch is the area between two
           consecutive `x` values. The value `Delta` corresponds to the distance between two
           consecutive values in x-direction. As we are only considering Cartesian grids, `Delta`
           is equal for all patches
- `Q`: Vector which contains the Control points
- `IP`: Coefficients matrix
"""
mutable struct CubicBSpline{x_type, Delta_type, Q_type, IP_type}
    x::x_type
    Delta::Delta_type
    Q::Q_type
    IP::IP_type
end

# Fill structure
@doc raw"""
    CubicBSpline(x::Vector, y::Vector; end_condition = "free", smoothing_factor = 0.0)

This function calculates the inputs for the structure [`CubicBSpline`](@ref).
The input values are:
- `x`: Vector that contains equally spaces values in x-direction
- `y`: Vector that contains values in y-direction
- `end_condition`: String that can either be `free` or `not-a-knot` and defines which
                   end condition should be considered. By default this is set to "free".
- `smoothing_factor`: Float64 ``\geq`` 0.0 which specifies the degree of smoothing of the `y`
                      values. By default this value is set to `0.0` that corresponds to no
                      smoothing.

First the data is sorted via [`sort_data`](@ref) to guarantee that the `x` values
are in ascending order.

The patch size `Delta` is calculated by subtracting the second and first `x` value. This can be done
because we only consider equally spaced `x` values.
(A patch is the area between two consecutive `x` values)

If a `smoothing_factor` > 0.0 is set, the function [`spline_smoothing`](@ref) calculates new `y` values
which guarantee a B-Spline with less curvature.

The coefficients matrix `IP` for linear B-splines is fixed to be
```math
\begin{aligned}
  IP = \begin{pmatrix}
    -1 & 3 & -3 & 1\\
    3 & -6 & 3 & 0\\
    -3 & 0 & 3 & 0\\
    1 & 4 & 1 & 0
  \end{pmatrix}
\end{aligned}
```

The "free" end condition requires the second and the second to last control points lie
between the first and the third control point and that the second to last control points are
between the third to last and the last control point. This procedure is only possible with at least two
values in `x` data. The system of linear equations to determine the control points have the following form:
```math
\begin{aligned}
    \underbrace{\begin{bmatrix}
            0 \\ P_1 \\ P_2 \\ \vdots \\ P_{n-1} \\ P_n\\ 0
        \end{bmatrix}}_{:= P^*_{\text{free}}}
        = \frac{1}{6}
    \underbrace{
        \begin{bmatrix}
            1 & -2 & 1 & 0 & ... & ... & 0 \\
            1      & 4 & 1 & 0 & ... & ... & 0\\
            0      & 1 & 4 & 1 & 0   &     &     \vdots\\
            \vdots &  0      & \ddots & \ddots & \ddots & 0 & \vdots\\
            \vdots &       & 0 & 1 & 4 & 1 & 0\\
            0 & ... & ... & 0 & 1 & 4 & 1\\
            0 & ... & ... & 0 & 1 & -2 & 1
        \end{bmatrix}
    }_{:= \Phi^*_{\text{free}}}
    \underbrace{\begin{bmatrix}
        Q_1 \\ Q_2 \\ Q_3 \\ \vdots \\ Q_n \\ Q_{n+1} \\ Q_{n+2}
    \end{bmatrix}}_{:= Q_{\text{free}}},
\end{aligned}
```
which is solved for ``Q_{\text{free}}``.

The "not-a-knot" end condition requires the continuity of the third derivative in the second
and second to last fit knot. This end condition is only possible with at least four values in `x` data.
The system of linear equations to determine the control points has the following form:
```math
\begin{aligned}
    \underbrace{\begin{bmatrix}
        0 \\ P_1 \\ P_2 \\ \vdots \\ P_{n-1} \\ P_n\\ 0
    \end{bmatrix}}_{:= P^*_{\text{not-a-knot}}}
    = \frac{1}{6}
    \underbrace{
        \begin{bmatrix}
            -1 & 4 & -6 & 4 & -1 & 0 &... &  0 \\
            1      & 4 & 1 & 0 & ... & ... & ... & 0\\
            0      & 1 & 4 & 1 & 0   &     &     &\vdots\\
            \vdots &  0      & \ddots & \ddots & \ddots & 0 & &\vdots\\
            \vdots &  & 0      & \ddots & \ddots & \ddots & 0 &\vdots\\
            \vdots &    &   & 0 & 1 & 4 & 1 & 0\\
            0 & ... & ...  & ... & 0 & 1 & 4 & 1\\
            0 & ... & 0 & -1 & 4 & -6 & 4 & -1
        \end{bmatrix}
    }_{:= \Phi^*_{\text{not-a-knot}}}
    \underbrace{\begin{bmatrix}
        Q_1 \\ Q_2 \\ Q_3 \\ \vdots \\ \vdots \\ Q_n \\ Q_{n+1} \\ Q_{n+2}
    \end{bmatrix}}_{:= Q_{\text{not-a-knot}}}.
\end{aligned}
```
which is solved for ``Q_{\text{not-a-knot}}``.

For both cases ``P_1,...,P_n = y_1,...,y_n``.

A reference for the calculations in this script can be found in Chapter 1 of
-  Quentin Agrapart & Alain Batailly (2020),
   Cubic and bicubic spline interpolation in Python.
   [hal-03017566v2](https://hal.archives-ouvertes.fr/hal-03017566v2)
"""
function CubicBSpline(x::Vector, y::Vector; end_condition = "free", smoothing_factor = 0.0)
    if length(x) < 2
        throw(ArgumentError("To perform cubic B-spline interpolation, we need an x vector
                             which contains at least 2 values."))
    end

    x, y = sort_data(x, y)

    Delta = x[2] - x[1]

    # Consider spline smoothing if required
    if smoothing_factor > 0
        y = spline_smoothing(smoothing_factor, Delta, y)
    end

    n = length(x)
    P = vcat(0, y, 0)
    IP = @SMatrix [-1 3 -3 1;
                   3 -6 3 0;
                   -3 0 3 0;
                   1 4 1 0]

    # Free end condition
    if end_condition == "free"
        du = vcat(-2, ones(n))
        dm = vcat(1, 4 * ones(n), 1)
        dl = vcat(ones(n), -2)

        Phi = Matrix(Tridiagonal(dl, dm, du))
        Phi[1, 3] = 1
        Phi[end, end - 2] = 1
        Phi_free = sparse(Phi)
        Q_free = 6 * (Phi_free \ P)

        CubicBSpline(x, Delta, Q_free, IP)

        # Not-a-knot end condition
    elseif end_condition == "not-a-knot"
        if length(x) < 4
            throw(ArgumentError("To perform cubic B-spline interpolation with not-a-knot
                                 end condition, we need an x vector which contains
                                 at least 4 values."))
        end

        du = vcat(4, ones(n))
        dm = vcat(-1, 4 * ones(n), -1)
        dl = vcat(ones(n), 4)

        Phi = Matrix(Tridiagonal(dl, dm, du))
        Phi[1, 3:5] = [-6 4 -1]
        Phi[end, (end - 4):(end - 2)] = [-1 4 -6]
        Phi_knot = sparse(Phi)

        Q_knot = 6 * (Phi_knot \ P)

        CubicBSpline(x, Delta, Q_knot, IP)

    else
        throw(ArgumentError("Only \"free\" and \"not-a-knot\" boundary conditions
                             are available!"))
    end
end

# Read from file
"""
    CubicBSpline(path::String; end_condition = "free", smoothing_factor = 0.0)

A function that reads in the `x` and `y` values for [`CubicBSpline`](@ref) from a .txt file.
The input values are:
- `path`: String of a path of the specific .txt file
- `end_condition`: String which can either be `free` or `not-a-knot` and defines which
                   end condition should be considered.
                   By default this is set to `free`
- `smoothing_factor`: Float64 ``\\geq`` 0.0 which specifies the degree of smoothing of the `y` values.
                      By default this value is set to `0.0` which corresponds to no smoothing.

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
function CubicBSpline(path::String; kwargs...)
    x, y = parse_txt_1D(path)

    CubicBSpline(x, y; kwargs...)
end

########################
### Lavery Spline ###
########################

# Helper container so that the JuMP `model` only needs constructed once
mutable struct LaverySplineModel
    model::Model
    b::Vector{VariableRef}
    abs_b::Vector{VariableRef}
    abs_E::Matrix{VariableRef}
    deltaZ::Vector{Float64} # data-dependent coefficients
end

# Constructor
function build_lavery_spline_model(len::Int, weight::Float64, integral_steps::Int)
    sumDomain = 1:len - 1
    bDomain   = 1:len

    # integration grid
    integralDomain = collect(range(-0.5, 0.5; length = integral_steps))

    # precompute t-dependent scalars
    aa = -1 .+ 6 .* integralDomain
    bb =  1 .+ 6 .* integralDomain
    cc = 12 .* integralDomain

    model = direct_model(
        optimizer_with_attributes(
            HiGHS.Optimizer,
            "presolve" => "on",
            "solver"   => "ipm",
        )
    )
    set_silent(model)

    # variables
    @variable(model, b[bDomain])
    @variable(model, abs_b[bDomain] >= 0)
    @variable(model, abs_E[i in sumDomain, k in 1:integral_steps] >= 0)

    # placeholder for data
    deltaZ = zeros(len - 1)

    inv_steps = 1 / integral_steps

    @objective(model, Min,
        sum(inv_steps * abs_E[i, k] for i in sumDomain, k in 1:integral_steps) +
        sum(weight * abs_b[i] for i in bDomain)
    )

    # constraints
    @constraint(model, [i in sumDomain, k in 1:integral_steps],
        abs_E[i, k] >= aa[k] * b[i] + bb[k] * b[i + 1] - cc[k] * deltaZ[i]
    )

    @constraint(model, [i in sumDomain, k in 1:integral_steps],
        abs_E[i, k] >= -aa[k] * b[i] - bb[k] * b[i + 1] + cc[k] * deltaZ[i]
    )

    @constraint(model, [i in bDomain], abs_b[i] >=  b[i])
    @constraint(model, [i in bDomain], abs_b[i] >= -b[i])

    return LaverySplineModel(model,
                             b, abs_b, abs_E,
                             deltaZ)
end

# Lavery Spline structure
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
mutable struct LaverySpline1D{x_type, z_type, b_type}
    x::x_type
    z::z_type
    b::b_type
    weight::Float64
    integral_steps::Int
end

# Fill structure
@doc raw"""
    LaverySpline1D(xData::AbstractVector, zData::AbstractVector;
                   weight::Float64 = 1e-4,
                   integral_steps::Int = 10)

This function calculates the inputs for the structure [`LaverySpline`](@ref).
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
function LaverySpline1D(xData::Vector, zData::Vector;
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
    xData, zData = sort_data(collect(xData), collect(zData))

    len = length(xData)

    # Build the JuMP model once to save time
    spline_model = build_lavery_spline_model(len, weight, integral_steps)

    for i in 1:len-1
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

# Read from file
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
function LaverySpline1D(path::String; kwargs...)
    x, z = parse_txt_1D(path)

    LaverySpline1D(x, z; kwargs...)
end