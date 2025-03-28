
####################################################
### Helper functions for 1D spline interpolation ###
####################################################

# Sorting the inputs so that the x values are ascending
"""
    sort_data(x::Vector,y::Vector)

Sorts the input vectors `x` and `y` so that `x` is in ascending order and the `y` values
accordingly to still correspond to the `x` values.
"""
function sort_data(x::Vector,y::Vector)

  original_data = hcat(x,y)
  sorted_data = original_data[sortperm(original_data[:,1]), :]

  x_sorted = sorted_data[:,1]
  y_sorted = sorted_data[:,2]

  return x_sorted,y_sorted
end

# Spline smoothing
@doc raw"""
    spline_smoothing(lambda::Number, Delta::Number, y::Vector)

The inputs to this function are:
- `lambda`: Smoothing factor which specifies the degree of the smoothing that should take place
- `Delta`: Step size of a patch (A patch is the area between two consecutive `x` values)
- `y`: Data values to be smoothed

The goal is to find a new interpolation values ``\hat{y}`` for ``y``, so that for given ``\lambda``,
the following equation is minimized:
```math
\begin{aligned}
		\text{PSS} = \sum_{i = 1}^{n} \left( y_i - \underbrace{S(t_i)}_{=\hat{y}_i} \right)^2
    + \lambda \int_{x_1}^{x_n} (S''(t))^2 dt,
	\end{aligned}
```
where ``S(t)`` is a cubic spline function.
``\hat{y}`` is determined as follows:
```math
\begin{aligned}
\hat{y} = (I+\lambda K)^{-1} y
\end{aligned}
```
where ``I`` is the ``n \times n`` identity matrix and ``K = \Delta_2^T W^{-1} \Delta_2`` with
```math
\begin{aligned}
\Delta_2 = \begin{pmatrix}
1/\Delta & -2/\Delta & 1/\Delta & ... & 0\\
0 & \ddots & \ddots & \ddots & 0\\
0 & ... & 1/\Delta & -2/\Delta & 1/\Delta
\end{pmatrix} \in \mathbb{R}^{(n-2) \times n}
\end{aligned}
```
and
```math
\begin{aligned}
W = \begin{pmatrix}
2/3 \Delta & 1/6 \Delta & 0 & ... & 0\\
1/6 \Delta & 2/3 \Delta & 1/6 \Delta & ... & 0\\
0 & \ddots & \ddots & \ddots & 0\\
0 & ... & 0 & 2/3 \Delta & 1/6 \Delta
\end{pmatrix} \in \mathbb{R}^{n \times n}
\end{aligned}
```

- Germán Rodríguez (2001),
  [Smoothing and non-parametric regression](https://data.princeton.edu/eco572/smoothing.pdf)
"""
function spline_smoothing(lambda::Number, Delta::Number, y::Vector)

  n = length(y)

  Delta_vec = repeat([Delta], n-2)

  Delta_2_ii   =  1 ./ Delta_vec
  Delta_2_iip1 = -2 ./ Delta_vec
  Delta_2_iip2 =  1 ./ Delta_vec

  Delta_2             =  zeros(n-2, n)
  Delta_2[:, 1:(n-2)] =  diagm(Delta_2_ii)
  Delta_2[:, 2:(n-1)] += diagm(Delta_2_iip1)
  Delta_2[:, 3: n   ] += diagm(Delta_2_iip2)

  W_im1i =  Delta_vec[1:end-1] ./ 6
  W_ii   = (2*Delta_vec) ./ 3
  W      = SymTridiagonal(W_ii, W_im1i)

  K = transpose(Delta_2) * inv(W) * Delta_2

  return inv(diagm(ones(n)) + lambda*K) * y
end

