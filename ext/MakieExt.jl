# Package extension for adding Makie-based features to TrixiBottomTopography.jl
module MakieExt

# Required for visualization code
using Makie: Makie

# Use all exported symbols to avoid having to rewrite all the visualization routines
using TrixiBottomTopography

# Import functions such that they can be extended with new methods
import TrixiBottomTopography: evaluate_two_dimensional_interpolant, plot_topography,
                              plot_topography_with_interpolation_knots

# To keep functionality we allow for the old function name `evaluate_bicubicspline_interpolant` to be used.
# From a front end perspective nothing changes; however, if one starts Julia with `--depwarn=yes`
# a warning would be thrown when one uses the deprecated function name.
function TrixiBottomTopography.evaluate_bicubicspline_interpolant(spline, x, y)
    Base.depwarn("evaluate_bicubicspline_interpolant is deprecated. Use evaluate_two_dimensional_interpolant instead.",
                 :evaluate_bicubicspline_interpolant)
    return evaluate_two_dimensional_interpolant(spline, x, y)
end

"""
    evaluate_two_dimensional_interpolant(spline, x, y)

Helper function to sample the bicubic spline function `spline`
at the interpolation nodes `x` (with size `n`) and `y` (with size `m`)
and return the values in an array `z` of size `m` by `n`.
"""
function evaluate_two_dimensional_interpolant(spline, x, y)

    # Get dimensions for solution matrix
    n = length(x)
    m = length(y)

    # Create empty solution matrix
    z = zeros(m, n)

    # Evaluate spline function at given x, y values
    for i in 1:n, j in 1:m
        z[j, i] = spline(x[i], y[j])
    end

    # Return the interpolated values
    return z
end

# Set the Makie default colormap to `:greenbrownterrain`.
default_Makie_colormap() = :greenbrownterrain

"""
    plot_topography(x, y;
                    xlabel = "", ylabel= "",
                    color, legend_position = :rt)

Plot function for 1D bottom topography data.
The interpolated values are provided in the vectors `x` and `y`.
Axis labels can be prescribed with `xlabel`, `ylabel`, and `zlabel`.
The default `color` for the interpolated bottom topography is light blue.
The `legend_position` controls the placement of the legend.
"""
function plot_topography(x, y; xlabel = "", ylabel = "", color = Makie.wong_colors(2)[5])
    Makie.lines(x, y, axis = (xlabel = xlabel, ylabel = ylabel), color = color)
end

"""
    plot_topography_with_interpolation_knots(x, y, x_knots, y_knots;
                                             xlabel = "", ylabel= "",
                                             [color,]
                                             legend_position = :rt)

Plot function for 1D bottom topography data together with the interpolation knots.
The interpolated values are provided in the vectors `x` and `y`.
The interpolation knots are provided in a similar fashion in vectors `x_knots` and
`y_knots`.
Axis labels can be prescribed with `xlabel`, `ylabel`, and `zlabel`.
The default `color` for the interpolated bottom topography is light blue.
The `legend_position` controls the placement of the legend.
"""
function plot_topography_with_interpolation_knots(x, y, x_knots, y_knots;
                                                  xlabel = "", ylabel = "",
                                                  color = Makie.wong_colors(2)[5],
                                                  legend_position = :rt)
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1], xlabel = xlabel, ylabel = ylabel)

    # Approximate topography with a line
    Makie.lines!(ax, x, y, color = color, label = "Bottom topography")
    # Interpolation knots as orange dots
    Makie.scatter!(ax, x_knots, y_knots, color = :orange, strokewidth = 1,
                   strokecolor = (:black, 0.5), markersize = 5, label = "Knots")
    # Create the legend
    Makie.axislegend(ax, position = legend_position)
    # Show the figure
    fig
end

"""
    plot_topography(x, y, z;
                    xlabel = "", ylabel= "", zlabel = "",
                    colormap,
                    azimuth_angle = 1.275*pi,
                    elevation_angle = pi/8)

Plot function for 2D bottom topography data.
The interpolated values are provided in the vectors `x`, `y` and an array `z`
where the shape of array `z` must be `size(y)` by `size(x)`.
Axis labels can be prescribed with `xlabel`, `ylabel`, and `zlabel`.
The default `colormap` is `:greenbrownterrain`.
Other possible colormaps for bottom topographies are `:darkterrain`, `:sandyterrain`,
or `:gist_earth`.
The left / right camera angle is controlled by `azimuth_angle` whereas the up / down
camera angle is controlled by `elevation_angle`. Both angle arguments must be given in radians.
"""
function plot_topography(x, y, z; xlabel = "", ylabel = "", zlabel = "",
                         colormap = default_Makie_colormap(),
                         azimuth_angle = 1.275 * pi, elevation_angle = pi / 8)
    # Generic plot function for 2D bottom topography data

    Makie.surface(x, y, z,
                  axis = (type = Makie.Axis3,
                          xlabel = xlabel, ylabel = ylabel, zlabel = zlabel,
                          azimuth = azimuth_angle, elevation = elevation_angle),
                  colormap = colormap)
end

"""
    plot_topography_with_interpolation_knots(x, y, z, x_knots, y_knots, z_knots;
                                             xlabel = "", ylabel= "", zlabel = "",
                                             colormap, legend_position = :rt,
                                             azimuth_angle = 1.275*pi,
                                             elevation_angle = pi/8)

Plot function for 2D bottom topography data together with the interpolation knots.
The interpolated values are provided in the vectors `x`, `y` and an array `z`
where the shape of array `z` must have `size(y)` by `size(x)`.
The interpolation knots are provided in a similar fashion in vectors `x_knots`,
`y_knots` and an array `z_knots` where the shape of `z_knots` must have `size(y_knots)` by `size(x_knots)`.
Axis labels can be prescribed with `xlabel`, `ylabel`, and `zlabel`.
The default `colormap` is `:greenbrownterrain`.
Other possible colormaps for bottom topographies are `:darkterrain`, `:sandyterrain`,
or `:gist_earth`.
The `legend_position` controls the placement of the legend.
The left / right camera angle is controlled by `azimuth_angle` whereas the up / down
camera angle is controlled by `elevation_angle`. Both angle arguments must be given in radians.
"""
function plot_topography_with_interpolation_knots(x, y, z, x_knots, y_knots, z_knots;
                                                  xlabel = "", ylabel = "", zlabel = "",
                                                  colormap = default_Makie_colormap(),
                                                  legend_position = :rt,
                                                  azimuth_angle = 1.275 * pi,
                                                  elevation_angle = pi / 8)
    fig = Makie.Figure()
    ax = Makie.Axis3(fig[1, 1], xlabel = xlabel, ylabel = ylabel, zlabel = zlabel,
                     azimuth = azimuth_angle, elevation = elevation_angle)

    # Approximate topography with a line
    Makie.surface!(ax, x, y, z, colormap = colormap, label = "Bottom topography")
    # Interpolation knots as orange dots
    Makie.scatter!(ax, x_knots, y_knots, z_knots, color = :orange, strokewidth = 1,
                   strokecolor = (:black, 0.5), markersize = 5, label = "Knots")
    # Create the legend
    Makie.axislegend(ax, position = legend_position)
    # Show the figure
    fig
end

end
