
#######################################
### Bilinear B Spline interpolation ###
#######################################

# Bilinear B Spline structure
"""
    BilinearBSpline(x, y, Delta, Q, IP)

Two dimensional bilinear B-spline structure that contains all important attributes to define
a B-Spline interpolation function.
These attributes are:
- `x`: Vector of values in x-direction
- `y`: Vector of values in y-direction
- `Delta`: Length of one side of a single patch in the given data set. A patch is the area between
           two consecutive `x` and `y` values. The value `Delta` corresponds to the distance
           between two consecutive values in x-direction. As we are only considering Cartesian
           grids, `Delta` is equal for all patches in x and y-direction
- `Q`: Matrix which contains the control points
- `IP`: Coefficients matrix
"""
mutable struct BilinearBSpline{x_type, y_type, Delta_type, Q_type, IP_type}
    x::x_type
    y::y_type
    Delta::Delta_type
    Q::Q_type
    IP::IP_type
end

# Fill structure
@doc raw"""
    BilinearBSpline(x::Vector, y::Vector, z::Matrix)

This function calculates the inputs for the structure [`BilinearBSpline`](@ref).
The input values are:
- `x`: Vector that contains equally spaced values in x-direction
- `y`: Vector that contains equally spaced values in y-direction where the spacing between the
       y-values has to be the same as the spacing between the x-values
- `z`: Matrix that contains the corresponding values in z-direction.
       Where the values are ordered in the following way:
```math
\begin{aligned}
   \begin{matrix}
        & & x_1 & x_2 & ... & x_n\\
        & & & & &\\
        y_1 & & z_{11} & z_{12} & ... & z_{1n}\\
        y_1 & & z_{21} & z_{22} & ... & z_{2n}\\
        \vdots & & \vdots & \vdots & \ddots & \vdots\\
        y_m & & z_{m1} & z_{m2} & ... & z_{mn}
   \end{matrix}
\end{aligned}
```
Bilinear B-spline interpolation is only possible if we have at least two values in `x`
and two values in `y` and the dimensions of vectors `x` and `y` correspond with the dimensions
of the matrix `z`.

First of all the data is sorted which is done by
[`sort_data`](@ref) to guarantee
that the `x` and `y` values are in ascending order with corresponding matrix `z`.

The patch size `Delta` is calculated by subtracting the second by the first `x` value. This can be
done because we only consider equal space between consecutive `x` and `y` values.
A patch is the area between two consecutive `x` and `y` values.

For bilinear B-spline interpolation, the control points `Q` correspond with the `z` values.

The coefficients matrix `IP` for bilinear B-splines is fixed to be
```math
\begin{aligned}
  \begin{pmatrix}
    -1 & 1\\
    1 & 0
  \end{pmatrix}
\end{aligned}
```

A reference for the calculations in this script can be found in Chapter 2 of
- Quentin Agrapart & Alain Batailly (2020),
  Cubic and bicubic spline interpolation in Python.
  [hal-03017566v2](https://hal.archives-ouvertes.fr/hal-03017566v2)
"""
function BilinearBSpline(x::Vector, y::Vector, z::Matrix)
    n = length(x)
    m = length(y)

    if (size(z, 2) != n || size(z, 1) != m)
        throw(DimensionMismatch("The dimensions of z do not coincide with x and y."))
    end

    if (n < 2 || m < 2)
        throw(ArgumentError("To perform bilinear B-spline interpolation, we need x and y
                             vectors which contain at least 2 values each."))
    end

    x, y, z = sort_data(x, y, z)

    Delta = x[2] - x[1]

    P = vcat(reshape(z', (m * n, 1)))
    IP = @SMatrix [-1 1;
                   1 0]

    Q = reshape(P, (n, m))

    BilinearBSpline(x, y, Delta, Q, IP)
end

# Read from file
"""
    BilinearBSpline(path::String)

A function which reads in the `x`, `y` and `z` values for
[`BilinearBSpline`](@ref) from a .txt
file. The input values are:
- `path`: String of a path of the specific .txt file

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
- Remaining lines: values for `z` where each value has its own line and is in the following order:
                   z_11, z_12, ... z_1n, z_21, ... z_2n, ..., z_m1, ..., z_mn

An example can be found [here](https://gist.githubusercontent.com/maxbertrand1996/7b1b943eac142d5bc836bb818fe83a5a/raw/74228e349e91fbfe1563479f99943b469f26ac62/Rhine_data_2D_10.txt)
"""
function BilinearBSpline(path::String)
    x, y, z = parse_txt_2D(path)

    BilinearBSpline(x, y, z)
end

######################################
### Bicubic B Spline interpolation ###
######################################

# Bicubic B-spline structure
"""
    BicubicBSpline(x, y, Delta, Q, IP)

Two dimensional cubic B-spline structure that contains all important attributes to define
a B-Spline interpolation function.
These attributes are:
- `x`: Vector of values in x-direction
- `y`: Vector of values in y-direction
- `Delta`: Length of one side of a single patch in the given data set. A patch is the area between
           two consecutive `x` and `y` values. The value `Delta` corresponds to the distance between
           two consecutive values in x-direction. As we are only considering Cartesian grids,
           `Delta` is equal for all patches in x and y-direction
- `Q`: Matrix which contains the control points
- `IP`: Coefficients matrix
"""
mutable struct BicubicBSpline{x_type, y_type, Delta_type, Q_type, IP_type}
    x::x_type
    y::y_type
    Delta::Delta_type
    Q::Q_type
    IP::IP_type
end

# Fill structure
@doc raw"""
    BicubicBSpline(x::Vector, y::Vector, z::Matrix; end_condition = "free", smoothing_factor = 0.0)

This function calculates the inputs for the structure [`BicubicBSpline`](@ref).
The input values are:
- `x`: Vector that contains equally spaced values in x-direction
- `y`: Vector that contains equally spaced values in y-direction where the spacing between the
       y-values has to be the same as the spacing between the x-values
- `z`: Matrix that contains the corresponding values in z-direction.
       Where the values are ordered in the following way:
```math
\begin{aligned}
    \begin{matrix}
        & & x_1 & x_2 & ... & x_n\\
        & & & & &\\
        y_1 & & z_{11} & z_{12} & ... & z_{1n}\\
        y_1 & & z_{21} & z_{22} & ... & z_{2n}\\
        \vdots & & \vdots & \vdots & \ddots & \vdots\\
        y_m & & z_{m1} & z_{m2} & ... & z_{mn}
    \end{matrix}
\end{aligned}
```
- `end_condition`: a string which can either be "free" or "not-a-knot" and defines which
                   end condition should be considered. By default this is set to "free".
- `smoothing_factor`: a Float64 ``\geq`` 0.0 which specifies the degree of smoothing of the `z`
                      values. By default this value is set to 0.0 which corresponds to no
                      smoothing.

Bicubic B-spline interpolation is only possible if the dimensions of vectors `x` and `y` correspond
with the dimensions of the matrix `z`.

First, the data is sorted via [`sort_data`](@ref) to guarantee
that the `x` and `y` values are in ascending order with corresponding matrix `z`.

The patch size `Delta` is calculated by subtracting the second by the first `x` value. This can be
done because we only consider equal space between consecutive `x` and `y` values.
A patch is the area between two consecutive `x` and `y` values.

If a `smoothing_factor` ``>`` 0.0 is set, the function [`calc_tps`](@ref)
calculates new values for `z` which guarantee a resulting parametric B-spline surface
with less curvature.

The coefficients matrix `IP` for bicubic B-splines is fixed to be
```math
\begin{aligned}
    \begin{pmatrix}
        -1 & 3 & -3 & 1\\
        3 & -6 & 3 & 0\\
        -3 & 0 & 3 & 0\\
        1 & 4 & 1 & 0
    \end{pmatrix}
\end{aligned}
```
To get the matrix of control points `Q` which is necessary to set up an interpolation function,
we need to define a matrix `Phi` which maps the control points to a vector `P`. This can be done
by solving the following system of linear equations for `Q`.
```math
\underbrace{
    \begin{bmatrix}
        z_{1,1} \\ z_{1,2} \\ \vdots \\ z_{1,n} \\ z_{2,1} \\ \vdots \\ z_{m,n} \\ 0 \\ \vdots \\ 0
    \end{bmatrix}
    }_{\text{:=P} \in \mathbb{R}^{(m+2)(n+2)\times 1}} = \frac{1}{36}
    \Phi \cdot
    \underbrace{\begin{bmatrix}
        Q_{1,1} \\ Q_{1,2} \\ \vdots \\ Q_{1,n+2} \\ Q_{2,1} \\ \vdots \\ Q_{m+2,n+2}
    \end{bmatrix}}_{\text{:= Q} \in \mathbb{R}^{(m+2) \times (n+2)}}
```
For the first `n` ``\cdot`` `m` lines, the matrix `Phi` is the same for the "free" end and the
"not-a-knot" end condition. These lines have to address the following condition:
```math
\begin{align*}
    z_{j,i} = \frac{1}{36} \Big( &Q_{j,i} + 4Q_{j+1,i} + Q_{j+2,i} + 4Q_{j,i+1} + 16Q_{j+1,i+1}\\
    &+ 4Q_{j+2,i+1} + Q_{j,i+2} + 4Q_{j+1,i+2} + Q_{j+2,i+2} \Big)
\end{align*}
```
for i = 1,...,n and j = 1,...,m.

The "free" end condition needs at least two values for the `x` and `y` vectors.
The "free" end condition has the following additional requirements for the control points which have
to be addressed by `Phi`:
- ``Q_{j,1} - 2Q_{j,2} + Q_{j,3} = 0`` for j = 2,...,m+1
- ``Q_{j,n} - 2Q_{j,n+1} + Q_{j,n+2} = 0`` for j = 2,...,m+1
- ``Q_{1,i} - 2Q_{2,i} + Q_{3,i} = 0`` for i = 2,...,n+1
- ``Q_{m,i} - 2Q_{m+1,i} + Q_{m+2,i} = 0`` for i = 2,...,n+1
- ``Q_{1,1} - 2Q_{2,2} + Q_{3,3} = 0``
- ``Q_{m+2,1} - 2Q_{m+1,2} + Q_{m,3} = 0``
- ``Q_{1,n+2} - 2Q_{2,n+1} + Q_{3,n} = 0``
- ``Q_{m,n} - 2Q_{m+1,n+1} + Q_{m+2,n+2} = 0``

The "not-a-knot" end condition needs at least four values for the `x` and `y` vectors.
- Continuity of the third `x` derivative between the leftmost and second leftmost patch
- Continuity of the third `x` derivative between the rightmost and second rightmost patch
- Continuity of the third `y` derivative between the patch at the top and the patch below
- Continuity of the third `y` derivative between the patch at the bottom and the patch above
- ``Q_{1,1} - Q_{1,2} - Q_{2,1} + Q_{2,2} = 0 ``
- ``Q_{m-1,1} + Q_{m,1} + Q_{m-1,2} - Q_{m,2} = 0``
- ``Q_{1,n-1} + Q_{2,n} + Q_{1,n-1} - Q_{2,n} = 0``
- ``Q_{m-1,n-1} - Q_{m,n-1} - Q_{m-1,n} + Q_{m,n} = 0``

A reference for the calculations in this script can be found in Chapter 2 of
- Quentin Agrapart & Alain Batailly (2020),
  Cubic and bicubic spline interpolation in Python.
  [hal-03017566v2](https://hal.archives-ouvertes.fr/hal-03017566v2)
"""
function BicubicBSpline(x::Vector, y::Vector, z::Matrix; end_condition = "free",
                        smoothing_factor = 0.0)
    n = length(x)
    m = length(y)

    if (size(z, 2) != n || size(z, 1) != m)
        throw(DimensionMismatch("The dimensions of z do not coincide with x and y."))
    end

    x, y, z = sort_data(x, y, z)

    # Consider spline smoothing if required
    if smoothing_factor > 0
        z = calc_tps(smoothing_factor, x, y, z)
    end

    Delta = x[2] - x[1]
    boundary_elmts = 4 + 2 * m + 2 * n
    inner_elmts = m * n
    P = vcat(reshape(z', (inner_elmts, 1)), zeros(boundary_elmts))

    IP = @SMatrix [-1 3 -3 1;
                   3 -6 3 0;
                   -3 0 3 0;
                   1 4 1 0]

    # Mapping matrix Phi
    Phi = spzeros((m + 2) * (n + 2), (m + 2) * (n + 2))

    # Fill inner control point matrix
    idx = 0
    for i in 1:inner_elmts
        Phi[i, idx + 1] = 1
        Phi[i, idx + 2] = 4
        Phi[i, idx + 3] = 1
        Phi[i, idx + (n + 2) + 1] = 4
        Phi[i, idx + (n + 2) + 2] = 16
        Phi[i, idx + (n + 2) + 3] = 4
        Phi[i, idx + 2 * (n + 2) + 1] = 1
        Phi[i, idx + 2 * (n + 2) + 2] = 4
        Phi[i, idx + 2 * (n + 2) + 3] = 1

        if (i % n) == 0
            idx += 3
        else
            idx += 1
        end
    end

    ########################
    ## Free end condition ##
    ########################
    if end_condition == "free"
        if (n < 2) || (m < 2)
            throw(ArgumentError("To perform bicubic B-spline interpolation with
                                 the free end condition, we need x and y vectors
                                 which contain at least 2 values each."))
        end

        # Q_{j,1} - 2Q_{j,2} + Q_{j,3} = 0
        idx = 0
        for i in (inner_elmts + 1):(inner_elmts + m)
            Phi[i, idx + (n + 2) + 1] = 1
            Phi[i, idx + (n + 2) + 2] = -2
            Phi[i, idx + (n + 2) + 3] = 1
            idx += (n + 2)
        end

        # Q_{j,n} - 2Q_{j,n+1} + Q_{j,n+2} = 0
        idx = 0
        for i in (inner_elmts + (m + 1)):(inner_elmts + (2 * m))
            Phi[i, idx + (n + 2) + (n)] = 1
            Phi[i, idx + (n + 2) + (n + 1)] = -2
            Phi[i, idx + (n + 2) + (n + 2)] = 1
            idx += (n + 2)
        end

        # Q_{1,i} - 2Q_{2,i} + Q_{3,i} = 0
        idx = 0
        for i in (inner_elmts + (2 * m) + 1):(inner_elmts + (2 * m) + n)
            Phi[i, idx + 2] = 1
            Phi[i, idx + (n + 2) + 2] = -2
            Phi[i, idx + 2 * (n + 2) + 2] = 1
            idx += 1
        end

        # Q_{m,i} - 2Q_{m+1,i} + Q_{m+2,i} = 0
        idx = (m - 1) * (n + 2)
        for i in (inner_elmts + (2 * m) + (n + 1)):(inner_elmts + (2 * m) + (2 * n))
            Phi[i, idx + 2] = 1
            Phi[i, idx + (n + 2) + 2] = -2
            Phi[i, idx + 2 * (n + 2) + 2] = 1
            idx += 1
        end

        i = inner_elmts + boundary_elmts - 3

        # Q_{1,1} - 2Q_{2,2} + Q_{3,3} = 0
        Phi[(i), 1] = 1
        Phi[(i), (n + 2) + 2] = -2
        Phi[(i), (2) * (n + 2) + 3] = 1

        # Q_{m+2,1} - 2Q_{m+1,2} + Q_{m,3} = 0
        Phi[(i + 1), (n + 2)] = 1
        Phi[(i + 1), (2) * (n + 2) - 1] = -2
        Phi[(i + 1), (3) * (n + 2) - 2] = 1

        # Q_{1,n+2} - 2Q_{2,n+1} + Q_{3,n} = 0
        Phi[(i + 2), (m - 1) * (n + 2) + 3] = 1
        Phi[(i + 2), (m) * (n + 2) + 2] = -2
        Phi[(i + 2), (m + 1) * (n + 2) + 1] = 1

        # Q_{m,n} - 2Q_{m+1,n+1} + Q_{m+2,n+2} = 0
        Phi[(i + 3), (m) * (n + 2) - 2] = 1
        Phi[(i + 3), (m + 1) * (n + 2) - 1] = -2
        Phi[(i + 3), (m + 2) * (n + 2)] = 1

        # For the sparse matrix `Phi` using the built-in `qr` function
        # is more numerically stable than the standard `LDLt` procedure
        # used by default in Julia for the `\` operator.
        Q_temp = 36 * (qr(Phi) \ P)
        Q = reshape(Q_temp, (n + 2, m + 2))

        BicubicBSpline(x, y, Delta, Q, IP)

        ###################################
        ## Not-a-knot boundary condition ##
        ###################################
    elseif end_condition == "not-a-knot"
        if (n < 4) || (m < 4)
            throw(ArgumentError("To perform bicubic B-spline interpolation with
                                 the not-a-knot end condition, we need x and y vectors
                                 which contain at least 4 values each."))
        end

        # Continuity of the third `x` derivative between the leftmost and second leftmost patch
        idx = 0
        for i in (inner_elmts + 1):(inner_elmts + m)
            Phi[i, idx + 1] = -1
            Phi[i, idx + 2] = 4
            Phi[i, idx + 3] = -6
            Phi[i, idx + 4] = 4
            Phi[i, idx + 5] = -1
            Phi[i, idx + (n + 2) + 1] = -4
            Phi[i, idx + (n + 2) + 2] = 16
            Phi[i, idx + (n + 2) + 3] = -24
            Phi[i, idx + (n + 2) + 4] = 16
            Phi[i, idx + (n + 2) + 5] = -4
            Phi[i, idx + 2 * (n + 2) + 1] = -1
            Phi[i, idx + 2 * (n + 2) + 2] = 4
            Phi[i, idx + 2 * (n + 2) + 3] = -6
            Phi[i, idx + 2 * (n + 2) + 4] = 4
            Phi[i, idx + 2 * (n + 2) + 5] = -1
            idx += (n + 2)
        end

        # Continuity of the third `x` derivative between the rightmost and second rightmost patch
        idx = (n + 2) + 1
        for i in (inner_elmts + (m + 1)):(inner_elmts + (2 * m))
            Phi[i, idx - 5] = -1
            Phi[i, idx - 4] = 4
            Phi[i, idx - 3] = -6
            Phi[i, idx - 2] = 4
            Phi[i, idx - 1] = -1
            Phi[i, idx + (n + 2) - 5] = -4
            Phi[i, idx + (n + 2) - 4] = 16
            Phi[i, idx + (n + 2) - 3] = -24
            Phi[i, idx + (n + 2) - 2] = 16
            Phi[i, idx + (n + 2) - 1] = -4
            Phi[i, idx + 2 * (n + 2) - 5] = -1
            Phi[i, idx + 2 * (n + 2) - 4] = 4
            Phi[i, idx + 2 * (n + 2) - 3] = -6
            Phi[i, idx + 2 * (n + 2) - 2] = 4
            Phi[i, idx + 2 * (n + 2) - 1] = -1
            idx += (n + 2)
        end

        # Continuity of the third `y` derivative between the patch at the top and the patch below
        idx = 0
        for i in (inner_elmts + (2 * m) + 1):(inner_elmts + (2 * m) + n)
            Phi[i, idx + 1] = -1
            Phi[i, idx + 2] = -4
            Phi[i, idx + 3] = -1
            Phi[i, idx + (n + 2) + 1] = 4
            Phi[i, idx + (n + 2) + 2] = 16
            Phi[i, idx + (n + 2) + 3] = 4
            Phi[i, idx + 2 * (n + 2) + 1] = -6
            Phi[i, idx + 2 * (n + 2) + 2] = -24
            Phi[i, idx + 2 * (n + 2) + 3] = -6
            Phi[i, idx + 3 * (n + 2) + 1] = 4
            Phi[i, idx + 3 * (n + 2) + 2] = 16
            Phi[i, idx + 3 * (n + 2) + 3] = 4
            Phi[i, idx + 4 * (n + 2) + 1] = -1
            Phi[i, idx + 4 * (n + 2) + 2] = -4
            Phi[i, idx + 4 * (n + 2) + 3] = -1
            idx += 1
        end

        # Continuity of the third `y` derivative between the patch at the bottom and the patch above
        idx = (m - 3) * (n + 2)
        for i in (inner_elmts + (2 * m) + (n + 1)):(inner_elmts + (2 * m) + (2 * n))
            Phi[i, idx + 1] = -1
            Phi[i, idx + 2] = -4
            Phi[i, idx + 3] = -1
            Phi[i, idx + (n + 2) + 1] = 4
            Phi[i, idx + (n + 2) + 2] = 16
            Phi[i, idx + (n + 2) + 3] = 4
            Phi[i, idx + 2 * (n + 2) + 1] = -6
            Phi[i, idx + 2 * (n + 2) + 2] = -24
            Phi[i, idx + 2 * (n + 2) + 3] = -6
            Phi[i, idx + 3 * (n + 2) + 1] = 4
            Phi[i, idx + 3 * (n + 2) + 2] = 16
            Phi[i, idx + 3 * (n + 2) + 3] = 4
            Phi[i, idx + 4 * (n + 2) + 1] = -1
            Phi[i, idx + 4 * (n + 2) + 2] = -4
            Phi[i, idx + 4 * (n + 2) + 3] = -1
            idx += 1
        end

        i = inner_elmts + boundary_elmts - 3

        # Q_{1,1} - Q_{1,2} - Q_{2,1} + Q_{2,2} = 0
        Phi[(i), 1] = 1
        Phi[(i), 2] = -1
        Phi[(i), (n + 2) + 1] = -1
        Phi[(i), (n + 2) + 2] = 1

        # Q_{m-1,1} + Q_{m,1} + Q_{m-1,2} - Q_{m,2} = 0
        Phi[(i + 1), (n + 2) - 1] = -1
        Phi[(i + 1), (n + 2)] = 1
        Phi[(i + 1), 2 * (n + 2) - 1] = 1
        Phi[(i + 1), 2 * (n + 2)] = -1

        # Q_{1,n-1} + Q_{2,n} + Q_{1,n-1} - Q_{2,n} = 0
        Phi[(i + 2), (m) * (n + 2) + 1] = -1
        Phi[(i + 2), (m) * (n + 2) + 2] = 1
        Phi[(i + 2), (m + 1) * (n + 2) + 1] = 1
        Phi[(i + 2), (m + 1) * (n + 2) + 2] = -1

        # Q_{m-1,n-1} - Q_{m,n-1} - Q_{m-1,n} + Q_{m,n} = 0
        Phi[(i + 3), (m + 1) * (n + 2) - 1] = 1
        Phi[(i + 3), (m + 1) * (n + 2)] = -1
        Phi[(i + 3), (m + 2) * (n + 2) - 1] = -1
        Phi[(i + 3), (m + 2) * (n + 2)] = 1

        # For the sparse matrix `Phi` using the built-in `qr` function
        # is more numerically stable than the standard `LDLt` procedure
        # used by default in Julia for the `\` operator.
        Q_temp = 36 * (qr(Phi) \ P)
        Q = reshape(Q_temp, (n + 2, m + 2))

        BicubicBSpline(x, y, Delta, Q, IP)
    else
        throw(ArgumentError("Only \"free\" and \"not-a-knot\" boundary
                             conditions are available!"))
    end
end

# Read from file
"""
    BicubicBSpline(path::String; end_condition = "free", smoothing_factor = 0.0)

A function which reads in the `x`, `y` and `z` values for
[`BilinearBSpline`](@ref)
from a .txt file. The input values are:
- `path`: String of a path of the specific .txt file
- `end_condition`: String which can either be "free" or "not-a-knot" and defines which
                   end condition should be considered.
                   By default this is set to "free"
- `smoothing_factor`: Float64 ``\\geq`` 0.0 which specifies the degree of smoothing of the `y` values.
                      By default this value is set to 0.0 which corresponds to no smoothing.

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

An example can be found [here](https://gist.githubusercontent.com/maxbertrand1996/7b1b943eac142d5bc836bb818fe83a5a/raw/74228e349e91fbfe1563479f99943b469f26ac62/Rhine_data_2D_10.txt)
"""
function BicubicBSpline(path::String; kwargs...)
    x, y, z = parse_txt_2D(path)

    BicubicBSpline(x, y, z; kwargs...)
end

# Lavery Spline structure
"""
    LaverySpline2D(x, y, z, bx, by, lambda, N)

Two dimensional Lavery spline structure.

The attributes are:
- `x`: Vector of data points in x-direction (knot points)
- `y`: Vector of data points in x-direction (knot points)
- `z`: Matrix of data values in the z-direction
- `bx`: Matrix of slope coefficients at each knot point
- `by`: Matrix of slope coefficients at each knot point
- `lambda`: Regularization weight for the slope coefficients (default: 1.0)
- `N`: Number of integration steps for computation (default: 5)

"""
mutable struct LaverySpline2D
    x::Vector
    y::Vector
    z::Array
    bx::Array
    by::Array
    lambda::Float64
    N::Int
end

# Fill structure
@doc raw"""
    LaverySpline2D(xData::Vector, yData::Vector, zData::Matrix;
                   lambda::Float64 = 1.0,
                   N::Int = 5)

This function calculates the inputs for the structure [`LaverySpline2D`](@ref).
The input values are:
- `xData`: Vector of x-coordinates of the data points (knots)
- `yData`: Vector of y-coordinates of the data points (knots)
- `zData`: Matrix of z-coordinates (function values) at the data points
- `lambda`: Regularization parameter for smoothness (default: 1.0)
- `N`: Number of discrete points for integration (default: 5)

First the data is sorted via [`sort_data`](@ref) to guarantee that the `x`, `y`, and `z`
values are in ascending order.

The Lavery spline is computed by solving an optimization problem.

The optimization problem is formulated as a linear program and solved using HiGHS.jl.

References:
- John E. Lavery (2001),
   Shape-preserving, multiscale interpolation by bi- and multivariate cubic splines.
   [DOI: 10.1016/S0167-8396(01)00034-6](https://doi.org/10.1016/S0167-8396(01)00034-6)
- Logan Eriksson and Oscar Jemsson (2024)
   Uni- and bivariate interpolation of multiscale data using cubic L1 splines.
   [DiVA1918338](https://www.diva-portal.org/smash/get/diva2:1918338/FULLTEXT01.pdf)
"""
function LaverySpline2D(xData::Vector, yData::Vector, zData::Matrix;
                        lambda::Float64 = 1.0,
                        N::Int = 5)
    I = length(xData)
    J = length(yData)

    if (size(zData, 1) != I || size(zData, 2) != J)
        throw(DimensionMismatch("The dimensions of z do not coincide with x and y."))
    end

    xData, yData, zData = sort_data(xData, yData, zData)

    model = direct_model(
        optimizer_with_attributes(
            HiGHS.Optimizer,
            "presolve" => "on",
            "solver"   => "ipm",
        )
    )
    set_silent(model)

    deltaX = [xData[i + 1] - xData[i] for i in 1:I-1]   # x-step length
    deltaY = [yData[j + 1] - yData[j] for j in 1:J-1]   # y-step length

    # Convenience values for inverse of the grid spacing
    inv_dx = 1 ./ deltaX
    inv_dy = 1 ./ deltaY

    # N is the number of integration points in each sub-piece
    domain = collect(range(0, 1, step = 1 / (2 * N)))

    @variable(model, bx[1:I, 1:J]) # dzdx(x_i,y_i)
    @variable(model, by[1:I, 1:J]) # dzdy(x_i,y_i)

    # Compute all the second derivatives in each of the patches
    d2zdx2_1 = [inv_dx[i]^2 * ((-6 + 12 * xTilde) * zData[i, j]
                + deltaX[i] * (-4 + 6 * xTilde) * bx[i, j]
                + (6 - 12 * xTilde) * zData[i + 1, j]
                + deltaX[i] * (-2 + 6 * xTilde) * bx[i + 1, j])
                for i in 1:I-1, j in 1:J-1, xTilde in domain]

    d2zdxdy_1= [inv_dx[i] * inv_dy[j] * (
                6 * yTilde * ((zData[i, j] + zData[i + 1, j + 1]) - (zData[i + 1, j] + zData[i, j + 1]))
                +deltaX[i] * yTilde * ((bx[i, j] + bx[i + 1, j]) - (bx[i, j + 1] + bx[i + 1, j + 1]))
                +deltaY[j] * ((by[i + 1, j] - by[i, j]) +
                2 * yTilde * ( (by[i, j] + by[i, j + 1]) - (by[i + 1, j] + by[i + 1, j + 1]))))
                for i in 1:I-1, j in 1:J-1, yTilde in domain]

    d2zdy2_1 = [inv_dy[j]^2 * ((-6 + 6 * xTilde + 6 * yTilde) * zData[i, j]
                + deltaX[i] * (-1 + xTilde) * bx[i, j]
                + deltaY[j] * (-3 + 2 * xTilde + 3 * yTilde) * by[i, j]
                + (-6 * xTilde + 6 * yTilde) * zData[i + 1, j]
                + deltaX[i] * xTilde * bx[i + 1, j]
                + deltaY[j] * (-1 - 2 * xTilde + 3 * yTilde) * by[i + 1, j]
                + (6 - 6 * xTilde - 6 * yTilde) * zData[i, j + 1]
                + deltaX[i] * (1 - xTilde) * bx[i, j + 1]
                + deltaY[j] * (-2 + 2 * xTilde + 3 * yTilde) * by[i, j + 1]
                + (6 * xTilde - 6 * yTilde) * zData[i + 1, j + 1]
                + deltaX[i] * (-xTilde) * bx[i + 1, j + 1]
                + deltaY[j] * (-2 * xTilde + 3 * yTilde) * by[i + 1, j + 1])
                for i in 1:I-1, j in 1:J-1, xTilde in domain, yTilde in domain]

    d2zdx2_2 = [inv_dy[j]^2 * ((-6 + 12 * yTilde) * zData[i + 1, j]
                + deltaY[j] * (-4 + 6 * yTilde) * by[i + 1, j]
                + (6 - 12 * yTilde) * zData[i + 1, j + 1]
                + deltaY[j] * (-2 + 6 * yTilde) * by[i + 1, j + 1])
                for i in 1:I-1, j in 1:J-1, yTilde in domain]

    d2zdxdy_2 = [inv_dy[j] * inv_dx[i] * (
                6 * (1 - xTilde) * ((zData[i + 1, j] + zData[i, j + 1]) - (zData[i + 1, j + 1] + zData[i, j]))
                + deltaY[j] * (1 - xTilde) * ((by[i + 1, j] + by[i + 1, j + 1]) - (by[i, j] + by[i, j + 1]) )
                + deltaX[i] * ((-bx[i + 1, j + 1] + bx[i + 1, j]) +
                2 * (1 - xTilde) * ((bx[i + 1, j] + bx[i, j]) - (-bx[i + 1, j + 1] - bx[i, j + 1]))))
                for i in 1:I-1, j in 1:J-1, xTilde in domain]

    d2zdy2_2 = [inv_dx[i]^2 * ((-6 + 6 * yTilde + 6 * (1 - xTilde)) * zData[i + 1, j]
                + deltaY[j] * (-1 + yTilde) * by[i + 1, j]
                + deltaX[i] * (-3 + 2 * yTilde + 3 * (1 - xTilde)) * (-bx[i + 1, j])
                + (-6 * yTilde + 6 * (1 - xTilde)) * zData[i + 1, j + 1]
                + deltaY[j] * yTilde * by[i + 1, j + 1]
                + deltaX[i] * (-1 - 2 * yTilde + 3 * (1 - xTilde)) * (-bx[i + 1, j + 1])
                + (6 - 6 * yTilde - 6 * (1 - xTilde)) * zData[i, j]
                + deltaY[j] * (1 - yTilde) * by[i, j]
                + deltaX[i] * (-2 + 2 * yTilde + 3 * (1 - xTilde)) * (-bx[i,j])
                + (6 * yTilde - 6 * (1 - xTilde)) * zData[i, j + 1]
                + deltaY[j] * (-yTilde) * by[i, j + 1]
                + deltaX[i] * (-2 * yTilde + 3 * (1 - xTilde)) * (-bx[i, j + 1]))
                for i in 1:I-1, j in 1:J-1, xTilde in domain, yTilde in domain]

    d2zdx2_3 = [inv_dx[i]^2 * ((-6 + 12 * (1 - xTilde)) * zData[i + 1, j + 1]
                + deltaX[i] * (-4 + 6 * (1 - xTilde)) * bx[i + 1, j + 1]
                + (6 - 12 * (1 - xTilde)) * zData[i, j + 1]
                + deltaX[i] * (-2 + 6 * (1 - xTilde)) * bx[i, j + 1])
                for i in 1:I-1, j in 1:J-1, xTilde in domain]

    d2zdxdy_3 = [inv_dx[i] * inv_dy[j] * (
                6 * (1 - yTilde) * ((zData[i + 1, j + 1] + zData[i, j]) - (zData[i, j + 1] + zData[i + 1, j]))
                + deltaX[i] * (1 - yTilde) * ((bx[i + 1, j + 1] + bx[i, j + 1]) - (bx[i + 1, j] + bx[i, j]))
                + deltaY[j] * ((by[i, j + 1] - by[i + 1, j + 1]) +
                2 * (1 - yTilde) * ((by[i + 1, j + 1] + by[i + 1, j]) - (by[i, j + 1] + by[i, j]))))
                for i in 1:I-1, j in 1:J-1, yTilde in domain]

    d2zdy2_3 = [inv_dy[j]^2 * ((-6 + 6 * (1 - xTilde) + 6 * (1 - yTilde)) * zData[i + 1, j + 1]
                + deltaX[i] * (-1 + (1 - xTilde)) * bx[i + 1, j + 1]
                + deltaY[j] * (-3 + 2 * (1 - xTilde) + 3 * (1 - yTilde)) * by[i + 1, j + 1]
                + (-6 * (1 - xTilde) + 6 * (1 - yTilde)) * zData[i, j + 1]
                + deltaX[i] * (1 - xTilde) * bx[i, j + 1]
                + deltaY[j] * (-1 - 2 * (1 - xTilde) + 3 * (1 - yTilde)) * by[i, j + 1]
                + (6 - 6 * (1 - xTilde) - 6 * (1 - yTilde)) * zData[i + 1, j]
                + deltaX[i] * (1 - (1 - xTilde )) * bx[i + 1, j]
                + deltaY[j] * (-2 + 2 * (1 - xTilde) + 3 * (1 - yTilde)) * by[i + 1, j]
                + (6 * (1 - xTilde) - 6 * (1 - yTilde)) * zData[i, j]
                + deltaX[i] * (-(1 - xTilde)) * bx[i, j]
                + deltaY[j] * (-2 * (1 - xTilde) + 3 * (1 - yTilde)) * by[i, j])
                for i in 1:I-1, j in 1:J-1, xTilde in domain, yTilde in domain]

    d2zdx2_4 = [inv_dy[j]^2 * ((-6 + 12 * (1 - yTilde)) * zData[i, j + 1]
                + deltaY[j] * (-4 + 6 * (1 - yTilde)) * (-bx[i, j + 1])
                + (6 - 12 * (1 - yTilde)) * zData[i, j]
                + deltaY[j] * (-2 + 6 * (1 - yTilde)) * (-bx[i,j]))
                for i in 1:I-1, j in 1:J-1, yTilde in domain]

    d2zdxdy_4 = [inv_dx[i] * inv_dy[j] * (
                6 * xTilde * ((zData[i, j + 1] + zData[i + 1, j]) - (zData[i, j] + zData[i + 1, j + 1]))
                + deltaY[j] * xTilde * ((-bx[i, j + 1] - bx[i, j]) - (-bx[i + 1, j + 1] - bx[i + 1, j]))
                + deltaX[i] * ((by[i, j] - by[i, j + 1])
                + 2 * xTilde * ((by[i, j + 1] + by[i + 1, j + 1]) - (by[i, j] + by[i + 1, j]))))
                for i in 1:I-1, j in 1:J-1, xTilde in domain]

    d2zdy2_4 = [inv_dx[i]^2 * ((-6 + 6 * (1 - yTilde) + 6 * xTilde) * zData[i, j + 1]
                + deltaY[j] * (-1 + (1 - yTilde)) * (-bx[i, j + 1])
                + deltaX[i] * (-3 + 2 * (1 - yTilde) + 3 * xTilde) * by[i, j + 1]
                + (-6 * (1 - yTilde) + 6 * xTilde) * zData[i, j]
                + deltaY[j] * (1 - yTilde) * (-bx[i, j])
                + deltaX[i] * (-1 -2 * (1 - yTilde) + 3 * xTilde) * by[i, j]
                + (6 - 6 * (1 - yTilde) - 6 * xTilde) * zData[i + 1, j + 1]
                + deltaY[j] * (1 - (1 - yTilde)) * (-bx[i + 1, j + 1])
                + deltaX[i] * (-2 + 2 * (1 - yTilde) + 3 * yTilde) * by[i + 1, j + 1]
                + (6 * (1 - yTilde) - 6 * xTilde) * zData[i + 1, j]
                + deltaY[j] * (-(1 - yTilde)) * (-bx[i + 1, j])
                + deltaX[i] * (-2 * (1 - yTilde) + 3 * xTilde) * by[i + 1,j])
                for i in 1:I-1, j in 1:J-1, xTilde in domain, yTilde in domain]

    # Define the variables
    @variable(model, abs_d2zdx2_1[1:I-1, 1:J-1, 1:(2*N)] >= 0)
    @variable(model, abs_d2zdxdy_1[1:I-1, 1:J-1, 1:(2*N)] >= 0)
    @variable(model, abs_d2zdy2_1[1:I-1, 1:J-1, 1:(2*N), 1:(2*N)] >= 0)
    @variable(model, abs_d2zdx2_2[1:I-1, 1:J-1, 1:(2*N)] >= 0)
    @variable(model, abs_d2zdxdy_2[1:I-1, 1:J-1, 1:(2*N)] >= 0)
    @variable(model, abs_d2zdy2_2[1:I-1, 1:J-1, 1:(2*N), 1:(2*N)] >= 0)
    @variable(model, abs_d2zdx2_3[1:I-1, 1:J-1, 1:(2*N)] >= 0)
    @variable(model, abs_d2zdxdy_3[1:I-1, 1:J-1, 1:(2*N)] >= 0)
    @variable(model, abs_d2zdy2_3[1:I-1, 1:J-1, 1:(2*N), 1:(2*N)] >= 0)
    @variable(model, abs_d2zdx2_4[1:I-1, 1:J-1, 1:(2*N)] >= 0)
    @variable(model, abs_d2zdxdy_4[1:I-1, 1:J-1, 1:(2*N)] >= 0)
    @variable(model, abs_d2zdy2_4[1:I-1, 1:J-1, 1:(2*N), 1:(2*N)] >= 0)
    @variable(model, abs_bx[1:I, 1:J] >= 0)
    @variable(model, abs_by[1:I, 1:J] >= 0)

    # Objective function

    # The sample size in each small square is 4*N^2 and lambda is a small number
    inv_N2 = 1 / (N^2)
    @objective(model, Min, sum(sum(inv_N2 * (
            sum(sum(abs_d2zdx2_1[i, j, k] + 2 * abs_d2zdxdy_1[i, j, l] + abs_d2zdy2_1[i, j, k, l]
                  + abs_d2zdx2_2[i, j, l] + 2 * abs_d2zdxdy_2[i, j, k] + abs_d2zdy2_2[i, j, k, l]
                  + abs_d2zdx2_3[i, j, k] + 2 * abs_d2zdxdy_3[i, j, l] + abs_d2zdy2_3[i, j, k, l]
                  + abs_d2zdx2_4[i, j, l] + 2 * abs_d2zdxdy_4[i, j, k] + abs_d2zdy2_4[i, j, k, l]
                for l in 1:k) for k in 1:N)
          + sum(sum(abs_d2zdx2_1[i, j, k] + 2 * abs_d2zdxdy_1[i, j, l] + abs_d2zdy2_1[i, j, k, l]
                  + abs_d2zdx2_2[i, j, l] + 2 * abs_d2zdxdy_2[i, j, k] + abs_d2zdy2_2[i, j, k, l]
                  + abs_d2zdx2_3[i, j, k] + 2 * abs_d2zdxdy_3[i, j, l] + abs_d2zdy2_3[i, j, k, l]
                  + abs_d2zdx2_4[i, j, l] + 2 * abs_d2zdxdy_4[i, j, k] + abs_d2zdy2_4[i, j, k, l]
                for l in N:2*N-k) for k in N+1:2*N)
          ) for j in 1:J-1) for i in 1:I-1)
          + lambda * sum(sum(abs_bx[i, j] + abs_by[i, j] for i in 1:I) for j in 1:J))

    # Constraints in first piece
    @constraint(model, pos_d2zdx2_1[i in 1:I-1, j in 1:J-1, k in 1:(2*N)], abs_d2zdx2_1[i, j, k] >=  d2zdx2_1[i, j, k])
    @constraint(model, neg_d2zdx2_1[i in 1:I-1, j in 1:J-1, k in 1:(2*N)], abs_d2zdx2_1[i, j, k] >= -d2zdx2_1[i, j, k])

    @constraint(model, pos_d2zdxdy_1[i in 1:I-1, j in 1:J-1, l in 1:(2*N)], abs_d2zdxdy_1[i, j, l] >=  d2zdxdy_1[i, j, l])
    @constraint(model, neg_d2zdxdy_1[i in 1:I-1, j in 1:J-1, l in 1:(2*N)], abs_d2zdxdy_1[i, j, l] >= -d2zdxdy_1[i, j, l])

    @constraint(model, pos_d2zdy2_1[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N)], abs_d2zdy2_1[i, j, k, l] >=  d2zdy2_1[i, j, k, l])
    @constraint(model, neg_d2zdy2_1[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N)], abs_d2zdy2_1[i, j, k, l] >= -d2zdy2_1[i, j, k, l])

    # Constraints in second piece
    @constraint(model, pos_d2zdx2_2[i in 1:I-1, j in 1:J-1, l in 1:(2*N)], abs_d2zdx2_2[i, j, l] >=  d2zdx2_2[i, j, l])
    @constraint(model, neg_d2zdx2_2[i in 1:I-1, j in 1:J-1, l in 1:(2*N)], abs_d2zdx2_2[i, j, l] >= -d2zdx2_2[i, j, l])

    @constraint(model, pos_d2zdxdy_2[i in 1:I-1, j in 1:J-1, k in 1:(2*N)], abs_d2zdxdy_2[i, j, k] >=  d2zdxdy_2[i, j, k])
    @constraint(model, neg_d2zdxdy_2[i in 1:I-1, j in 1:J-1, k in 1:(2*N)], abs_d2zdxdy_2[i, j, k] >= -d2zdxdy_2[i, j, k])

    @constraint(model, pos_d2zdy2_2[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N)], abs_d2zdy2_2[i, j, k, l] >=  d2zdy2_2[i, j, k, l])
    @constraint(model, neg_d2zdy2_2[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N)], abs_d2zdy2_2[i, j, k, l] >= -d2zdy2_2[i, j, k, l])

    # Constraints in third piece
    @constraint(model, pos_d2zdx2_3[i in 1:I-1, j in 1:J-1, k in 1:(2*N)], abs_d2zdx2_3[i, j, k] >=  d2zdx2_3[i, j, k])
    @constraint(model, neg_d2zdx2_3[i in 1:I-1, j in 1:J-1, k in 1:(2*N)], abs_d2zdx2_3[i, j, k] >= -d2zdx2_3[i, j, k])

    @constraint(model, pos_d2zdxdy_3[i in 1:I-1, j in 1:J-1, k in 1:(2*N)], abs_d2zdxdy_3[i, j, k] >=  d2zdxdy_3[i, j, k])
    @constraint(model, neg_d2zdxdy_3[i in 1:I-1, j in 1:J-1, k in 1:(2*N)], abs_d2zdxdy_3[i, j, k] >= -d2zdxdy_3[i, j, k])

    @constraint(model, pos_d2zdy2_3[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N)], abs_d2zdy2_3[i, j, k, l] >=  d2zdy2_3[i, j, k, l])
    @constraint(model, neg_d2zdy2_3[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N)], abs_d2zdy2_3[i, j, k, l] >= -d2zdy2_3[i, j, k, l])

    # Constraints in fourth piece
    @constraint(model, pos_d2zdx2_4[i in 1:I-1, j in 1:J-1, l in 1:(2*N)], abs_d2zdx2_4[i, j, l] >=  d2zdx2_4[i, j, l])
    @constraint(model, neg_d2zdx2_4[i in 1:I-1, j in 1:J-1, l in 1:(2*N)], abs_d2zdx2_4[i, j, l] >= -d2zdx2_4[i, j, l])

    @constraint(model, pos_d2zdxdy_4[i in 1:I-1, j in 1:J-1, k in 1:(2*N)], abs_d2zdxdy_4[i, j, k] >=  d2zdxdy_4[i, j, k])
    @constraint(model, neg_d2zdxdy_4[i in 1:I-1, j in 1:J-1, k in 1:(2*N)], abs_d2zdxdy_4[i, j, k] >= -d2zdxdy_4[i, j, k])

    @constraint(model, pos_d2zdy2_4[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N)], abs_d2zdy2_4[i, j, k, l] >=  d2zdy2_4[i, j, k, l])
    @constraint(model, neg_d2zdy2_4[i in 1:I-1, j in 1:J-1, k in 1:(2*N), l in 1:(2*N)], abs_d2zdy2_4[i, j, k, l] >= -d2zdy2_4[i, j, k, l])

    # Remaining constraints on the coefficients
    @constraint(model, pos_bx[i in 1:I, j in 1:J], abs_bx[i, j] >=  bx[i, j])
    @constraint(model, neg_bx[i in 1:I, j in 1:J], abs_bx[i, j] >= -bx[i, j])

    @constraint(model, pos_by[i in 1:I, j in 1:J], abs_by[i, j] >=  by[i, j])
    @constraint(model, neg_by[i in 1:I, j in 1:J], abs_by[i, j] >= -by[i, j])

    # Solve the optimization problem
    optimize!(model)

    bx = Array(value.(bx))
    by = Array(value.(by))

    return LaverySpline2D(xData, yData, zData, bx, by, lambda, N)
end

# Read from file
"""
    LaverySpline2D(path::String; lambda::Float64 = 1.0, N::Int = 5)

A function which reads in the `x`, `y` and `z` values for [`LaverySpline2D`](@ref) from a .txt file.
The input values are:
- `path`: String of a path of the specific .txt file
- `lambda`: Regularization parameter for smoothness (default: 1.0)
- `N`: Number of discrete points for integration (default: 5)

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
    x, y, z = parse_txt_2D(path)

    # The Lavery spline constructor expects the `z` data in a slightly different order
    z = Matrix(transpose(z))
    LaverySpline2D(x, y, z; kwargs...)
end