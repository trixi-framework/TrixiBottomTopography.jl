
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
function spline_interpolation(lavery::LaverySpline2D, x::Number, y::Number)
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
