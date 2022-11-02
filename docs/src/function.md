# B-spline interpolation function

With the defined structures, we can use these to set up the B-spline interpolation functions. The function to do so uses function overloading and therefore has the same name for all cases. The cases are distinguished by the given structure.

1D:
```julia
julia> spline_func(x) = spline_interpolation(spline_struct, x)
```
2D:
```julia
julia> spline_func(x, y) = spline_interpolation(spline_struct, x, y)
```

For further information, see:
- [`spline_interpolation(linear)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.spline_interpolation-Tuple{TrixiBottomTopography.LinearBSpline,%20Any}) 
- [`spline_interpolation(cubic)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.spline_interpolation-Tuple{TrixiBottomTopography.CubicBSpline,%20Any})
- [`spline_interpolation(bilinear)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.spline_interpolation-Tuple{TrixiBottomTopography.BilinearBSpline,%20Any,%20Any})
- [`spline_interpolation(bicubic)`](https://maxbertrand1996.github.io/TrixiBottomTopography.jl/dev/reference/#TrixiBottomTopography.spline_interpolation-Tuple{TrixiBottomTopography.BicubicBSpline,%20Any,%20Any})
