
##############################################
### Two dimensional B-spline interpolation ###
##############################################

# Bilinaer B-spline interpolation
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
-  Quentin Agrapart & Alain Batailly (2020)
   Cubic and bicubic spline interpolation in Python.
   [hal-03017566v2](https://hal.archives-ouvertes.fr/hal-03017566v2)
"""
function spline_interpolation(b_spline::BilinearBSpline, x::Number, y::Number)

  x_vec = b_spline.x
  y_vec = b_spline.y
  Delta = b_spline.Delta
  Q     = b_spline.Q
  IP    = b_spline.IP

  i = max(1, min(searchsortedlast(x_vec, x), length(x_vec) - 1))
  j = max(1, min(searchsortedlast(y_vec, y), length(y_vec) - 1))

  my = (x - x_vec[i])/Delta
  ny = (y - y_vec[j])/Delta

  Q_temp = [Q[i, j:(j+1)] Q[(i+1), j:(j+1)]]

  c = [ny, 1]' * IP * Q_temp * IP' * [my, 1]

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
  Q     = b_spline.Q
  IP    = b_spline.IP

  i = max(1, min(searchsortedlast(x_vec, x), length(x_vec) - 1))
  j = max(1, min(searchsortedlast(y_vec, y), length(y_vec) - 1))

  my = (x - x_vec[i])/Delta
  ny = (y - y_vec[j])/Delta

  Q_temp = [Q[i, j:(j+3)] Q[(i+1), j:(j+3)] Q[(i+2), j:(j+3)] Q[(i+3), j:(j+3)]]

  c = 1/36 * [ny^3, ny^2, ny, 1]' * IP * Q_temp * IP' * [my^3, my^2, my, 1]

  return c
end
