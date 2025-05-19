
##############################################
### One dimensional B-spline interpolation ###
##############################################

# Linear B-spline interpolation
@doc raw"""
    spline_interpolation(b_spline::LinearBSpline, x::Number)

The inputs are the [`LinearBSpline`](@ref) object and a variable `x` at which the spline
will be evaluated.

The parameter `i` indicates the patch in which the variable `x` is located.
This parameter is also used to get the correct control points from `Q`.
A patch is the area between two consecutive `b_spline.x` values.

`kappa` is an interim variable which maps `x` to the interval ``[0,1]``
for further calculations.

To evaluate the spline at `x`, we have to calculate the following:
```math
\begin{aligned}
c_{i,1}(\kappa_i(x)) =
    \begin{bmatrix}
        \kappa_i(x)\\ 1
    \end{bmatrix}^T
    \begin{bmatrix}
        -1 & 1\\1 & 0
    \end{bmatrix}
    \begin{bmatrix}
        Q_i\\Q_{i+1}
    \end{bmatrix}
\end{aligned}
```

A reference for the calculations in this script can be found in Chapter 1 of
- Quentin Agrapart & Alain Batailly (2020),
  Cubic and bicubic spline interpolation in Python.
  [hal-03017566v2](https://hal.archives-ouvertes.fr/hal-03017566v2)
"""
function spline_interpolation(b_spline::LinearBSpline, x::Number)
    x_vec = b_spline.x
    Delta = b_spline.Delta
    Q = b_spline.Q
    IP = b_spline.IP

    i = max(1, min(searchsortedlast(x_vec, x), length(x_vec) - 1))

    kappa_i = (x - x_vec[i]) / Delta

    # Allocation free version of the naive implementation
    # c_i1 = [kappa_i, 1]' * IP * Q[i:(i + 1)]

    # Note, IP has already been constructed as an `SMatrix`
    # in the `LinearBSpline` constructor
    kappa_vec = @SVector [kappa_i, 1]
    Q_slice = @SVector [Q[i], Q[i+1]]

    c_i1 = kappa_vec' * (IP * Q_slice)

    return c_i1
end

# Cubic B-spline interpolation
@doc raw"""
    spline_interpolation(b_spline::CubicBSpline, x::Number)

The inputs are the [`CubicBSpline`](@ref) object and a variable `x` at which the spline
will be evaluated.

The parameter `i` indicates the patch in which the variable `x` is located.
This parameter is also used to get the correct control points from `Q`.
A patch is the area between two consecutive `b_spline.x` values.

`kappa` is  an interim variable which maps `t` to the interval ``[0,1]``
for further calculations.

To evaluate the spline at `x`, we have to calculate the following:
```math
\begin{aligned}
c_{i,3}\left(\kappa_i(x) \right) = \frac{1}{6}
    \begin{bmatrix}
        \kappa_i(x)^3\\ \kappa_i(x)^2\\ \kappa_i(x) \\1
    \end{bmatrix}^T
    \underbrace{\begin{bmatrix}
        -1 & 3 & -3 & 1\\
        3 & -6 & 3 & 0\\
        -3 & 0 & 3 & 0\\
        1 & 4 & 1 & 0
    \end{bmatrix}}_{\text{IP}}
    \begin{bmatrix}
        Q_{i,\text{free}}\\ Q_{i+1,\text{free}}\\ Q_{i+2,\text{free}}\\ Q_{i+3,\text{free}}
    \end{bmatrix}
\end{aligned}
```

A reference for the calculations in this script can be found in Chapter 1 of
- Quentin Agrapart & Alain Batailly (2020),
  Cubic and bicubic spline interpolation in Python.
  [hal-03017566v2](https://hal.archives-ouvertes.fr/hal-03017566v2)
"""
function spline_interpolation(b_spline::CubicBSpline, x::Number)
    x_vec = b_spline.x
    Delta = b_spline.Delta
    Q = b_spline.Q
    IP = b_spline.IP

    i = max(1, min(searchsortedlast(x_vec, x), length(x_vec) - 1))

    kappa_i = (x - x_vec[i]) / Delta

    # Allocation free version of the naive implementation
    # c_i3 = 1 / 6 * [kappa_i^3, kappa_i^2, kappa_i, 1]' * IP * Q[i:(i + 3)]

    # Note, IP has already been constructed as an `SMatrix`
    # in the `CubicBSpline` constructor
    kappa_vec = @SVector [kappa_i^3, kappa_i^2, kappa_i, 1]
    Q_slice = @SVector [Q[i], Q[i+1], Q[i+2], Q[i+3]]

    c_i3 = (1 / 6) * (kappa_vec' * (IP * Q_slice))

    return c_i3
end
