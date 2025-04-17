# Package extension for adding Makie-based features to TrixiBottomTopography.jl
module TrixiBottomTopographyMakieExt

# Required for visualization code
using Makie: Makie

# Use all exported symbols to avoid having to rewrite all the visualization routines
using TrixiBottomTopography

# Import functions such that they can be extended with new methods
import TrixiBottomTopography: evaluate_bicubicspline_interpolant, plot_topography, plot_topography_with_interpolation_knots


# Helper function to sample the bicubic spline function `spline`
# at the nodes `x` and `y` and return the values in `z`
function evaluate_bicubicspline_interpolant(spline, x, y)

  # Get dimensions for solution matrix
  n = length(x)
  m = length(y)

  # Create empty solution matrix
  z = zeros(n, m)

  # Evaluate spline function at given x, y values
  for i in 1:n, j in 1:m
    z[j,i] = spline(x[i], y[j])
  end

  # Return the interpolated values
  return z
end

# We set the Makie default colormap to `:greenbrownterrain`.
# Other colormaps that may be of use are `:darkterrain`, `:sandyterrain` or `:gist_earth`
# Notes, Plots.jl default is `:inferno`.
default_Makie_colormap() = :greenbrownterrain

# Generic plot function for 1D bottom topography data
function plot_topography(x, y; xlabel = "", ylabel= "", color = Makie.wong_colors(2)[5])
  Makie.lines(x, y, axis=(xlabel = xlabel, ylabel = ylabel), color = color)
end

# Plot function for 1D bottom topography data together with the interpolation knots
function plot_topography_with_interpolation_knots(x, y, x_knots, y_knots;
                                                  xlabel = "", ylabel= "",
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

# Generic plot function for 2D bottom topography data
function plot_topography(x, y, z; xlabel = "", ylabel= "", zlabel = "",
                         colormap = default_Makie_colormap())
  Makie.surface(x, y, z, axis=(type = Makie.Axis3,
                xlabel = xlabel, ylabel = ylabel, zlabel = zlabel),
                colormap = colormap)
end

# Plot function for 2D bottom topography data together with the interpolation knots
function plot_topography_with_interpolation_knots(x, y, z, x_knots, y_knots, z_knots;
                                                  xlabel = "", ylabel= "", zlabel = "",
                                                  colormap = default_Makie_colormap(),
                                                  legend_position = :rt)
  fig = Makie.Figure()
  ax = Makie.Axis3(fig[1, 1], xlabel = xlabel, ylabel = ylabel, zlabel = zlabel)

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
