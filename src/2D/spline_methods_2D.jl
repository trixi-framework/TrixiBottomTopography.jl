
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
    slope = dy_j / dx_i

    # Depending in which Sibson element the `(x,y)` point is present the variables may change roles.
    if (y <= yData[j] + (x - xData[i]) * slope && y <= yData[j+1] - (x - xData[i]) * slope)
        # Sibson element 1 (default)
        xt = xTilde
        yt = yTilde
        dx = dx_i
        dy = dy_j
        z_ij = zData[i, j]
        zx_ij = bx[i, j]
        zy_ij = by[i, j]
        z_ip1j = zData[i + 1, j]
        zx_ip1j = bx[i + 1, j]
        zy_ip1j = by[i + 1, j]
        z_ijp1 = zData[i, j + 1]
        zx_ijp1 = bx[i, j + 1]
        zy_ijp1 = by[i, j + 1]
        z_ip1jp1 = zData[i + 1, j + 1]
        zx_ip1jp1 = bx[i + 1, j + 1]
        zy_ip1jp1 = by[i + 1, j + 1]
    elseif (y <= yData[j] + (x - xData[i]) * slope && y >= yData[j+1] - (x - xData[i]) * slope)
        # Sibson element 2; Eq. (2.27)
        xt = yTilde
        yt = 1 - xTilde
        dx = dy_j
        dy = dx_i
        z_ij = zData[i + 1, j]
        zx_ij = by[i + 1, j]
        zy_ij = -bx[i + 1, j]
        z_ip1j = zData[i + 1, j + 1]
        zx_ip1j = by[i + 1, j + 1]
        zy_ip1j = -bx[i + 1, j + 1]
        z_ijp1 = zData[i, j]
        zx_ijp1 = by[i, j]
        zy_ijp1 = -bx[i, j]
        z_ip1jp1 = zData[i, j + 1]
        zx_ip1jp1 = by[i, j + 1]
        zy_ip1jp1 = -bx[i, j + 1]
    elseif (y >= yData[j] + (x - xData[i]) * slope && y >= yData[j+1] - (x - xData[i]) * slope)
        # Sibson element 3; Eq. (2.28)
        xt = 1 - xTilde
        yt = 1 - yTilde
        dx = dx_i
        dy = dy_j
        z_ij = zData[i + 1, j + 1]
        zx_ij = -bx[i + 1, j + 1]
        zy_ij = -by[i + 1, j + 1]
        z_ip1j = zData[i, j + 1]
        zx_ip1j = -bx[i, j + 1]
        zy_ip1j = -by[i, j + 1]
        z_ijp1 = zData[i + 1, j]
        zx_ijp1 = -bx[i + 1, j]
        zy_ijp1 = -by[i + 1, j]
        z_ip1jp1 = zData[i, j]
        zx_ip1jp1 = -bx[i, j]
        zy_ip1jp1 = -by[i, j]
    elseif (y >= yData[j] + (x - xData[i]) * slope && y <= yData[j+1] - (x - xData[i]) * slope)
        # Sibson element 4; Eq. (2.29)
        xt = 1 - yTilde
        yt = xTilde
        dx = dy_j
        dy = dx_i
        z_ij = zData[i, j + 1]
        zx_ij = -by[i, j + 1]
        zy_ij = bx[i, j + 1]
        z_ip1j = zData[i, j]
        zx_ip1j = -by[i, j]
        zy_ip1j = bx[i, j]
        z_ijp1 = zData[i + 1, j + 1]
        zx_ijp1 = -by[i + 1, j + 1]
        zy_ijp1 = bx[i + 1, j + 1]
        z_ip1jp1 = zData[i + 1, j]
        zx_ip1jp1 = -by[i + 1, j]
        zy_ip1jp1 = bx[i + 1, j]
    end

    # Evaluation of the spline on Sibson element 1. This is Eq. (2.26) in Eriksson and Jemsson
    z_val = ((1 - 3 * xt^2 + 2 * xt^3 - 3 * yt^2 + 3 * xt * yt^2 + yt^3) * z_ij
            + dx * (xt - 2 * xt^2 + xt^3 - 0.5 * yt^2 + 0.5 * xt * yt^2) * zx_ij
            + dy * (yt - xt * yt - 1.5 * yt^2 + xt * yt^2 + 0.5 * yt^3) * zy_ij
            + (3 * xt^2 - 2 * xt^3 - 3 * xt * yt^2 + yt^3) * z_ip1j
            + dx * (-xt^2 + xt^3 + 0.5 * xt * yt^2) * zx_ip1j
            + dy * (xt * yt - 0.5 * yt^2 - xt * yt^2 + 0.5 * yt^3) * zy_ip1j
            + (3 * yt^2 - 3 * xt * yt^2 - yt^3) * z_ijp1
            + dx * (0.5 * yt^2 - 0.5 * xt * yt^2) * zx_ijp1
            + dy * (-yt^2 + xt * yt^2 + 0.5 * yt^3) * zy_ijp1
            + (3 * xt * yt^2 - yt^3) * z_ip1jp1
            + dx * (-0.5 * xt * yt^2) * zx_ip1jp1
            + dy * (-xt * yt^2 + 0.5 * yt^3) * zy_ip1jp1)

    return z_val
end
