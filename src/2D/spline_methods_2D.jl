
##############################################
### Two dimensional B-spline interpolation ###
##############################################

# Bilinear B-spline interpolation
@doc raw"""
    spline_interpolation(b_spline::BilinearBSpline, x::Number, y::Number)

The inputs are the [`BilinearBSpline`](@ref) object and the variable `x` and `y` at which the spline
will be evaluated.

The parameters `i` and `j` indicate the patch in which `(x,y)` is located.
This information is also used to get the correct control points from `Q`.
A patch is the area between two consecutive `b_spline.x` and `b_spline.y` values.

`my` is  an interim variable that maps `x` to the interval ``[0,1]``
for further calculations. `ny` does the same for `y`.

To evaluate the spline at `(x,y)`, we have to calculate the following:
```math
\begin{aligned}
c_{i,j,1}(\mu_i(x),\nu_j(y)) =
    \begin{bmatrix} \nu_j(y)\\ 1 \end{bmatrix}^T
    \underbrace{\begin{bmatrix} -1 & 1\\ 1 & 0 \end{bmatrix}}_{\text{IP}}
    \begin{bmatrix} Q_{i,j} & Q_{i,j+1}\\ Q_{i+1,j} & Q_{i+1,j+1} \end{bmatrix}
    \underbrace{\begin{bmatrix} -1 & 1\\ 1 & 0 \end{bmatrix}}_{\text{IP}^T}
    \begin{bmatrix} \mu_i(x) \\ 1\end{bmatrix}
\end{aligned}
```

A reference for the calculations in this script can be found in Chapter 2 of
- Quentin Agrapart & Alain Batailly (2020)
  Cubic and bicubic spline interpolation in Python.
  [hal-03017566v2](https://hal.archives-ouvertes.fr/hal-03017566v2)
"""
function spline_interpolation(b_spline::BilinearBSpline, x::Number, y::Number)
    x_vec = b_spline.x
    y_vec = b_spline.y
    Delta = b_spline.Delta
    Q = b_spline.Q
    IP = b_spline.IP

    i = max(1, min(searchsortedlast(x_vec, x), length(x_vec) - 1))
    j = max(1, min(searchsortedlast(y_vec, y), length(y_vec) - 1))

    my = (x - x_vec[i]) / Delta
    ny = (y - y_vec[j]) / Delta

    # Allocation free version of the naive implementation
    # Q_temp = [Q[i, j:(j + 1)] Q[(i + 1), j:(j + 1)]]
    # c = [ny, 1]' * IP * Q_temp * IP' * [my, 1]

    # Note, IP has already been constructed as an `SMatrix`
    # in the `BilinearBSpline` constructor
    ny_vec = @SVector [ny, 1]
    my_vec = @SVector [my, 1]
    Q_temp = @SMatrix [Q[i, j] Q[i + 1, j];
                       Q[i, j + 1] Q[i + 1, j + 1]]

    c = ny_vec' * IP * Q_temp * IP' * my_vec

    return c
end

# Bicubic B-spline interpolation
@doc raw"""
    spline_interpolation(b_spline::BicubicBSpline, x::Number, y::Number)

The inputs are the [`BicubicBSpline`](@ref) object and the variable `x` and `y` at which the spline
will be evaluated.

The parameters `i` and `j` indicate the patch in which `(x,y)` is located.
This information is also used to get the correct control points from `Q`.
A patch is the area between two consecutive `b_spline.x` and `b_spline.y` values.

`my` is  an interim variable which maps `x` to the interval ``[0,1]``
for further calculations. `ny` does the same for `y`.

To evaluate the spline at `(x,y)`, we have to calculate the following:
```math
\begin{aligned}
c_{i,j,3}(\mu_i(x),\nu_j(y)) = \frac{1}{36}
    \begin{bmatrix} \nu_j^3(y) \\ \nu_j^2(y) \\ \nu_j(y) \\ 1 \end{bmatrix}^T
    \underbrace{\begin{bmatrix}
        -1 & 3 & -3 & 1\\
        3 & -6 & 3 & 0\\
        -3 & 0 & 3 & 0\\
        1 & 4 & 1 & 0
    \end{bmatrix}}_{\text{IP}}
    \begin{bmatrix}
        Q_{i,j} & Q_{i+1,j} & Q_{i+2,j} & Q_{i+3,j}\\
        Q_{i,j+1} & Q_{i+1,j+1} & Q_{i+2,j+1} & Q_{i+3,j+1}\\
        Q_{i,j+2} & Q_{i+1,j+2} & Q_{i+2,j+2} & Q_{i+3,j+2}\\
        Q_{i,j+3} & Q_{i+1,j+3} & Q_{i+2,j+3} & Q_{i+3,j+3}
    \end{bmatrix}
    \underbrace{\begin{bmatrix}
        -1 & 3 & -3 & 1\\
        3 & -6 & 0 & 4\\
        -3 & 3 & 3 & 1\\
        1 & 0 & 0 & 0
    \end{bmatrix}}_{\text{IP}}
    \begin{bmatrix} \mu_i^3(x) \\ \mu_i^2(x) \\ \mu_i(x) \\ 1 \end{bmatrix}
\end{aligned}
```

A reference for the calculations in this script can be found in Chapter 2 of
-  Quentin Agrapart & Alain Batailly (2020),
   Cubic and bicubic spline interpolation in Python.
   [hal-03017566v2](https://hal.archives-ouvertes.fr/hal-03017566v2)
"""
function spline_interpolation(b_spline::BicubicBSpline, x::Number, y::Number)
    x_vec = b_spline.x
    y_vec = b_spline.y
    Delta = b_spline.Delta
    Q = b_spline.Q
    IP = b_spline.IP

    i = max(1, min(searchsortedlast(x_vec, x), length(x_vec) - 1))
    j = max(1, min(searchsortedlast(y_vec, y), length(y_vec) - 1))

    my = (x - x_vec[i]) / Delta
    ny = (y - y_vec[j]) / Delta

    # Allocation free version of the naive implementation
    # Q_temp = [Q[i, j:(j + 3)] Q[(i + 1), j:(j + 3)] Q[(i + 2), j:(j + 3)] Q[(i + 3), j:(j + 3)]]
    # c = 1 / 36 * [ny^3, ny^2, ny, 1]' * IP * Q_temp * IP' * [my^3, my^2, my, 1]

    # Note, IP has already been constructed as an `SMatrix`
    # in the `BicubicBSpline` constructor
    ny_vec = @SVector [ny^3, ny^2, ny, 1]
    my_vec = @SVector [my^3, my^2, my, 1]
    Q_temp = @SMatrix [Q[i, j] Q[i + 1, j] Q[i + 2, j] Q[i + 3, j];
                       Q[i, j + 1] Q[i + 1, j + 1] Q[i + 2, j + 1] Q[i + 3, j + 1];
                       Q[i, j + 2] Q[i + 1, j + 2] Q[i + 2, j + 2] Q[i + 3, j + 2];
                       Q[i, j + 3] Q[i + 1, j + 3] Q[i + 2, j + 3] Q[i + 3, j + 3]]

    c = (1 / 36) * (ny_vec' * IP * Q_temp * IP' * my_vec)

    return c
end

# Lavery spline interpolation
@doc raw"""
    spline_interpolation(lavery_spline::LaverySpline2D, x::Number, y::Number)

Evaluates the Bicubic Lavery spline at a single point `(x,y)`.
The spline on a rectangular patch is written in terms of four (triangular) Sibson elements.
Function form of the spline changes depending on which Sibson element contains the point `(x,y)`.
See Section 2.5 of the provided reference for details.

Reference:
- Logan Eriksson and Oscar Jemsson (2024)
   Uni- and bivariate interpolation of multiscale data using cubic L1 splines.
   [DiVA1918338](https://www.diva-portal.org/smash/get/diva2:1918338/FULLTEXT01.pdf)
"""
function spline_interpolation(lavery_spline::LaverySpline2D, x::Number, y::Number)
    xData = lavery_spline.x
    yData = lavery_spline.y
    zData = lavery_spline.z
    bx = lavery_spline.bx
    by = lavery_spline.by

    # Find the patch containing (x, y)
    i = max(1, min(searchsortedlast(xData, x), length(xData) - 1))
    j = max(1, min(searchsortedlast(yData, y), length(yData) - 1))

    # Width of patch
    dx_i = xData[i + 1] - xData[i]
    dy_j = yData[j + 1] - yData[j]

    # Local coordinate within patch
    xTilde = (x - xData[i]) / dx_i # Transforms x in [x_i,x_i+1] to xTilde in [0,1]
    yTilde = (y - yData[j]) / dy_j # Transforms y in [y_i,y_i+1] to yTilde in [0,1]

    z_val = 0
    slope = dy_j/dx_i

    # Depending on the location within the patch, use a different evaluation version
    if (y <= yData[j] + (x - xData[i]) * slope && y <= yData[j+1] - (x - xData[i]) * slope)
        z_val = ((1 - 3 * xTilde^2 + 2 * xTilde^3 - 3 * yTilde^2 + 3 * xTilde * yTilde^2 + yTilde^3) * zData[i, j]
                + dx_i * (xTilde - 2 * xTilde^2 + xTilde^3 - 0.5 * yTilde^2 + 0.5 * xTilde * yTilde^2) * bx[i, j]
                + dy_j * (yTilde - xTilde * yTilde - 1.5 * yTilde^2 + xTilde * yTilde^2 + 0.5 * yTilde^3) * by[i, j]
                + (3 * xTilde^2 - 2 * xTilde^3 - 3 * xTilde * yTilde^2 + yTilde^3) * zData[i + 1, j]
                + dx_i * (-xTilde^2 + xTilde^3 + 0.5 * xTilde * yTilde^2) * bx[i + 1, j]
                + dy_j * (xTilde * yTilde - 0.5 * yTilde^2 - xTilde * yTilde^2 + 0.5 * yTilde^3) * by[i + 1, j]
                + (3 * yTilde^2 - 3 * xTilde * yTilde^2 - yTilde^3) * zData[i, j + 1]
                + dx_i * (0.5 * yTilde^2 - 0.5 * xTilde * yTilde^2) * bx[i, j + 1]
                + dy_j * (-yTilde^2 + xTilde * yTilde^2 + 0.5 * yTilde^3) * by[i, j + 1]
                + (3 * xTilde * yTilde^2 - yTilde^3) * zData[i + 1, j + 1]
                + dx_i * (-0.5 * xTilde * yTilde^2) * bx[i + 1,j + 1]
                + dy_j * (-xTilde * yTilde^2 + 0.5 * yTilde^3) * by[i + 1, j + 1])
    elseif (y <= yData[j] + (x - xData[i]) * slope && y >= yData[j+1] - (x - xData[i]) * slope)
        z_val = ((1 - 3 * yTilde^2 + 2 * yTilde^3 - 3 * (1 - xTilde)^2 + 3 * yTilde * (1 - xTilde)^2 + (1 - xTilde )^3) * zData[i + 1, j]
                + dy_j * (yTilde - 2 * yTilde^2 + yTilde^3 - 0.5 * (1 - xTilde)^2 + 0.5 * yTilde * (1 - xTilde)^2) * by[i + 1, j]
                + dx_i * ((1 - xTilde) - yTilde * (1 - xTilde) - 1.5 * (1 - xTilde)^2 + yTilde * (1 - xTilde)^2 + 0.5 * (1 - xTilde)^3) * (-bx[i + 1, j])
                + (3 * yTilde^2 - 2 * yTilde^3 - 3 * yTilde * (1 - xTilde)^2 + (1 - xTilde)^3) * zData[i + 1, j + 1]
                + dy_j * (-yTilde^2 + yTilde^3 + 0.5 * yTilde * (1 - xTilde)^2) * by[i + 1, j + 1]
                + dx_i * (yTilde * (1 - xTilde) - 0.5 * (1 - xTilde)^2 - yTilde * (1 - xTilde)^2 + 0.5 * (1 - xTilde)^3) * (-bx[i + 1, j + 1])
                + (3 * (1 - xTilde)^2 - 3 * yTilde * (1 - xTilde)^2 - (1 - xTilde)^3) * zData[i, j]
                + dy_j * (0.5 * (1 - xTilde)^2 - 0.5 * yTilde * (1 - xTilde)^2) * by[i, j]
                + dx_i * (-(1 - xTilde)^2 + yTilde * (1 - xTilde)^2 + 0.5 * (1 - xTilde)^3) * (-bx[i, j])
                + (3 * yTilde * (1 - xTilde)^2 - (1 - xTilde)^3) * zData[i, j + 1]
                + dy_j * (-0.5 * yTilde * (1 - xTilde)^2) * by[i, j + 1]
                + dx_i * (-yTilde * (1 - xTilde)^2 + 0.5 * (1 - xTilde)^3) * (-bx[i, j + 1]))
    elseif (y >= yData[j] + (x - xData[i]) * slope && y >= yData[j+1] - (x - xData[i]) * slope)
        z_val = ((1 - 3 * (1 - xTilde)^2 + 2 * (1 - xTilde)^3 - 3 * (1 - yTilde)^2 + 3 * (1 - xTilde) * (1 - yTilde)^2 + (1 - yTilde)^3) * zData[i + 1, j + 1]
                + dy_j * ((1 - xTilde) - 2 * (1 - xTilde)^2 + (1 - xTilde)^3 - 0.5 * (1 - yTilde)^2 + 0.5 * (1 - xTilde) * (1 - yTilde)^2) * (-bx[i + 1, j + 1])
                + dx_i * ((1 - yTilde) - (1 - xTilde) * (1 - yTilde) - 1.5 * (1 - yTilde)^2 + (1 - xTilde) * (1 - yTilde)^2 + 0.5 * (1 - yTilde)^3) * (-by[i + 1, j + 1])
                + (3 * (1 - xTilde)^2 - 2 * (1 - xTilde)^3 - 3 * (1 - xTilde) * (1 - yTilde)^2 + (1 - yTilde)^3) * zData[i, j + 1]
                + dy_j * (-(1 - xTilde)^2 + (1 - xTilde)^3 + 0.5 * (1 - xTilde) * (1 - yTilde)^2) * (-bx[i, j + 1])
                + dx_i * ((1 - xTilde) * (1 - yTilde) - 0.5*(1 - yTilde)^2 - (1 - xTilde) * (1 - yTilde)^2 + 0.5 * (1 - yTilde)^3) * (-by[i, j + 1])
                + (3 * (1 - yTilde)^2 - 3 * (1 - xTilde) * (1 - yTilde)^2 - (1 - yTilde)^3) * zData[i + 1, j]
                + dy_j * (0.5 * (1 - yTilde)^2 - 0.5 * (1 - xTilde) * (1 - yTilde)^2) * (-bx[i + 1, j])
                + dx_i * (-(1 - yTilde)^2 + (1 - xTilde) * (1 - yTilde)^2 + 0.5 * (1 - yTilde)^3) * (-by[i+1,j])
                + (3 * (1 - xTilde) * (1 - yTilde)^2 - (1 - yTilde)^3) * zData[i, j]
                + dy_j * (-0.5 * (1 - xTilde) * (1 - yTilde)^2) * (-bx[i, j])
                + dx_i * (-(1 - xTilde) * (1 - yTilde)^2 + 0.5 * (1 - yTilde)^3) * (-by[i, j]))
    elseif (y >= yData[j] + (x - xData[i]) * slope && y <= yData[j+1] - (x - xData[i]) * slope)
        z_val = ((1 - 3 * (1 - yTilde)^2 + 2 * (1 - yTilde)^3 - 3 * xTilde^2 + 3 * (1 - yTilde) * xTilde^2 + xTilde^3) * zData[i, j + 1]
                + dy_j * ((1 - yTilde) - 2 * (1 - yTilde)^2 + (1 - yTilde)^3 - 0.5 * xTilde^2 + 0.5 * (1 - yTilde) * xTilde^2) * (-by[i, j + 1])
                + dx_i * (xTilde - (1 - yTilde) * xTilde - 1.5 * xTilde^2 + (1 - yTilde) * xTilde^2 + 0.5 * xTilde^3) * bx[i, j + 1]
                + (3 * (1 - yTilde)^2 - 2 * (1 - yTilde)^3 - 3 * (1 - yTilde) * xTilde^2 + xTilde^3) * zData[i, j]
                + dy_j * (-(1 - yTilde)^2 + (1 - yTilde)^3 + 0.5 * (1 - yTilde) * xTilde^2) * (-by[i, j])
                + dx_i * ((1 - yTilde) * xTilde - 0.5 * xTilde^2 - (1 - yTilde) * xTilde^2 + 0.5 * xTilde^3) * (bx[i, j])
                + (3 * xTilde^2 - 3 * (1 - yTilde) * xTilde^2 - xTilde^3) * zData[i + 1, j + 1]
                + dy_j * (0.5 * xTilde^2 - 0.5 * (1 - yTilde) * xTilde^2) * (-by[i + 1, j + 1])
                + dx_i * (-xTilde^2 + (1 - yTilde) * xTilde^2 + 0.5 * xTilde^3) * bx[i + 1, j + 1]
                + (3 * (1 - yTilde) * xTilde^2 - xTilde^3) * zData[i + 1, j]
                + dy_j * (-0.5 * (1 - yTilde) * xTilde^2) * (-by[i + 1, j])
                + dx_i * (-(1 - yTilde) * xTilde^2 + 0.5 * xTilde^3) * bx[i + 1, j])
    end

    return z_val
end
