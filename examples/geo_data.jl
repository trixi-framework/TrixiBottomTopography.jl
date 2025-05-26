############################################################################
# Script which uses GeophysicalModelGenerator to load real topography data #
# into TrixiBottomTopography.                                              #
############################################################################

#import Pkg
#Pkg.activate(@__DIR__)
#Pkg.instantiate()

#  Pkg.add("GeophysicalModelGenerator")
#  Pkg.add("GMT")
#  Pkg.add("CSV")
#  Pkg.add("DataFrames")
#  Pkg.add("Downloads")
#  Pkg.add("JuliaFormatter")
#  Pkg.add("DelimitedFiles")
#  Pkg.add("PyPlot")
#  Pkg.add("CairoMakie")
#  Pkg.add("Trixi")
#  Pkg.add("OrdinaryDiffEq")

using GMT
using GeophysicalModelGenerator
using CSV
using DataFrames
using TrixiBottomTopography
using Downloads: download
using JuliaFormatter
using DelimitedFiles
using CairoMakie
using OrdinaryDiffEq
using Trixi

###################################
# Specify the limits of the topography data based on the coordinates
# The following longitudianl and transversal coordinates are choosen
# to cover the same area as the other dataset in the examples folder
lon_min = 6.963880
lon_max = 6.978499
lat_min = 50.947861
lat_max = 50.957095
limits = [lon_min, lon_max, lat_min, lat_max]

lon_mean = (lon_max + lon_min) / 2
lat_mean = (lat_min + lat_max) / 2

###################################
#Loading the topography data
#important: if you have a small area you might use hiher resolution: 
# Note: 
#you can chose between different resolutions of the topography data.
# | Dataset                 |   Resolution |   Description                                               |
# |:----------------        | ------------ | ----------------------------------------------------------- |
# | "@earth\\_relief\\_01s" |	1 arc sec 	 | SRTM tiles (14297 tiles, land only, 60S-60N) [NASA/USGS]    |
# | "@earth\\_relief\\_03s"	|   3 arc sec	 | SRTM tiles (14297 tiles, land only, 60S-60N) [NASA/USGS]    |
# | "@earth\\_relief\\_15s"	|  15 arc sec	 | SRTM15+ [David Sandwell, SIO/UCSD]                          |
# | "@earth\\_relief\\_30s"	|  30 arc sec	 | SRTM30+ [Becker et al., 2009, SIO/UCSD]                     |
# | "@earth\\_relief\\_01m"	|   1 arc min	 | ETOPO1 Ice surface [NEIC/NOAA]                              |
# | "@earth\\_relief\\_02m"	|   2 arc min	 | ETOPO2v2 Ice surface [NEIC/NOAA]                            |
# | "@earth\\_relief\\_03m"	|   3 arc min	 | ETOPO1 after Gaussian spherical filtering (5.6 km fullwidth)|
# | "@earth\\_relief\\_04m"	|   4 arc min	 | ETOPO1 after Gaussian spherical filtering (7.5 km fullwidth)|
# | "@earth\\_relief\\_05m"	|   5 arc min	 | ETOPO1 after Gaussian spherical filtering (9 km fullwidth)  |
# | "@earth\\_relief\\_06m"	|   6 arc min	 | ETOPO1 after Gaussia30n spherical filtering (10 km fullwidth) |
# | "@earth\\_relief\\_10m"	|  10 arc min	 | ETOPO1 after Gaussian spherical filtering (18 km fullwidth) |
# | "@earth\\_relief\\_15m"	|  20 arc min	 | ETOPO1 after Gaussian spherical filtering (28 km fullwidth) |
# | "@earth\\_relief\\_20m"	|  20 arc min	 | ETOPO1 after Gaussian spherical filtering (37 km fullwidth) |
# | "@earth\\_relief\\_30m"	|  30 arc min	 | ETOPO1 after Gaussian spherical filtering (55 km fullwidth) |
# | "@earth\\_relief\\_60m"	|  60 arc min	 | ETOPO1 after Gaussian spherical filtering (111 km fullwidth)|

Topo = import_topo(lon = [lon_min, lon_max], lat = [lat_min, lat_max],
                   file = "@earth_relief_01s") # here we load the topography data

p=ProjectionPoint(Lon = lon_mean, Lat = lat_mean) # to use cartasian coordiantes we choose a projection point here 

Topo_Cart = convert2CartData(Topo, p) # here we use the projection point to convert the topography data to cartesian coordinates
###################################
#now we have to define the gridpoints for the cartesian coordinates

low_x = -0.5
high_x = 0.499
gridsize_x = 0.001
low_y = -0.5
high_y = 0.499
gridzize_y = 0.001

values_x = collect(low_x:gridsize_x:high_x)
values_y = collect(low_y:gridzize_y:high_y)

###################################
# the gridpoints have to fullfill the condition: Int(sqrt(length)) in TrixiBottomTopography
# so check if we have a right choice for our gridpoints:
function safe_computation(values_x, values_y)
    try
        Int(sqrt(length(values_x) * length(values_y)))
        return Int(sqrt(length(values_x) * length(values_y))), false
    catch e
        if e isa InexactError
            println("There is an InexactError. the gridsize has to be adjusted")
            return nothing, nothing, true
        else
            rethrow(e)
        end
    end
end

safe_computation(values_x, values_y)

###################################

Topo_Cart_orth = CartData(xyz_grid(low_x:gridsize_x:high_x, low_y:gridzize_y:high_y, 0)) # create a grid

Topo_Cart_orth = project_CartData(Topo_Cart_orth, Topo, p) # project the topography data on the cartasian grid
###################################
# write the topography data to a file compatible with TrixiBottomTopography

df_x = DataFrame(Topo_Cart_orth.x.val[:, :, 1], :auto);

df_y = DataFrame(Topo_Cart_orth.y.val[:, :, 1], :auto);

df_z = DataFrame(Topo_Cart_orth.z.val[:, :, 1], :auto);

# combine x, y, z into a DataFrame, warning: you have to scale the values to meters
df_xyz = DataFrame(x = convert.(Float64, vec(Topo_Cart_orth.x.val[:, :, 1])) .* 1000,
                   y = convert.(Float64, vec(Topo_Cart_orth.y.val[:, :, 1])) .* 1000,
                   z = convert.(Float64, vec(Topo_Cart_orth.z.val[:, :, 1])) .* 1000)

data_dir = joinpath(@__DIR__, "data")

# Write directly to the data directory
output_file = joinpath(data_dir, "geo.xyz")
open(output_file, "w") do file
    for row in eachrow(df_xyz)
        rounded_row = [round(value, digits = 5) for value in row]
        println(file, join(rounded_row, " "))
    end
end

path_src_file = joinpath(data_dir, "geo.xyz")

path_out_file_1d_x = joinpath(data_dir, "rhine_data_1d_20_x_geo.txt")
path_out_file_1d_y = joinpath(data_dir, "rhine_data_1d_20_y_geo.txt")
path_out_file_2d = joinpath(data_dir, "rhine_data_2d_20_geo.txt")

# Convert data
convert_dgm_1d(path_src_file, path_out_file_1d_x; excerpt = 20, section = 100);
convert_dgm_1d(path_src_file, path_out_file_1d_y; excerpt = 20, direction = "y",
               section = 100)
convert_dgm_2d(path_src_file, path_out_file_2d; excerpt = 20)

#################
#now redoo the steps of the turoial and try if we can do b-spline interpolation and simulations with this data
#1D:

data = joinpath(data_dir, "rhine_data_1d_20_x_geo.txt")

# Define B-spline structure
spline_struct = CubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 999)

spline_func(x) = spline_interpolation(spline_struct, x)

# Define interpolation points
n = 100
x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))

# Get interpolated values
y_int_pts = spline_func.(x_int_pts)

plot_topography(x_int_pts, y_int_pts; xlabel = "x[m]", ylabel = "z[m]")

# Get the original interpolation knots
x_knots = spline_struct.x
y_knots = spline_func.(x_knots)

plot_topography_with_interpolation_knots(x_int_pts, y_int_pts, x_knots, y_knots;
                                         xlabel = "x[m]", ylabel = "z[m]")

equations = ShallowWaterEquations1D(gravity_constant = 1.0, H0 = 57.0)

# Defining initial condition for the dam break problem
function initial_condition_dam_break(x, t, equations::ShallowWaterEquations1D)
    inicenter = SVector(0.0)
    x_norm = x[1] - inicenter[1]
    r = abs(x_norm)

    # Calculate primitive variables
    H = r < 50 ? 70.0 : 57.0
    v = 0.0
    b = spline_func(x[1])

    return prim2cons(SVector(H, v, b), equations)
end

# Setting initial condition
initial_condition = initial_condition_dam_break

# Setting the boundary to be a reflective wall
boundary_condition = boundary_condition_slip_wall

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
solver = DGSEM(polydeg = 3, surface_flux = (flux_hll, flux_nonconservative_fjordholm_etal),
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the TreeMesh and setup a periodic mesh

coordinates_min = spline_struct.x[1]
coordinates_max = spline_struct.x[end]
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 3,
                n_cells_max = 10_000,
                periodicity = false)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solvers

tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# run the simulation

# define equidistant nodes in time for visualization of an animation
visnodes = range(tspan[1], tspan[2], length = 90)

# use a Runge-Kutta method with error-based time step size control
sol = solve(ode, RDPK3SpFSAL49(), abstol = 1.0e-8, reltol = 1.0e-8,
            saveat = visnodes)

# Create animation of the solution
j = Observable(1)
time = Observable(0.0)

pd_list = [PlotData1D(sol.u[i], semi) for i in 1:length(sol.t)]
f = Figure()
title_string = lift(t -> "time t = $(round(t, digits=3))", time)
ax = Axis(f[1, 1], xlabel = "x[m]", ylabel = "z[m]", title = title_string)

height = lift(i -> pd_list[i].data[:, 1], j)
bottom = lift(i -> pd_list[i].data[:, 3], j)
CairoMakie.lines!(ax, pd_list[1].x, height)
CairoMakie.lines!(ax, pd_list[1].x, bottom)
ylims!(ax, 30, 90)

record(f, "animation.gif", 1:length(pd_list)) do tt
    j[] = tt
    time[] = sol.t[tt]
end

####################
#now redoo the steps of the turoial and try if we can do b-spline interpolation and simulations with this data
#2D:

data = joinpath(data_dir, "rhine_data_2d_20_geo.txt")

spline_struct = BicubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 999)

# Define B-spline interpolation function
spline_func(x, y) = spline_interpolation(spline_struct, x, y)

# Define interpolation points
n = 100
x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))
y_int_pts = Vector(LinRange(spline_struct.y[1], spline_struct.y[end], n))

# Get interpolated matrix
z_int_pts = evaluate_bicubicspline_interpolant(spline_func, x_int_pts, y_int_pts)

plot_topography(x_int_pts, y_int_pts, z_int_pts;
                xlabel = "x\n [m]",
                ylabel = "y\n [m]",
                zlabel = "z\n [m]",
                azimuth_angle = 54 * pi / 180,
                elevation_angle = 27 * pi / 180)

# Get the original interpolation knots
x_knots = spline_struct.x
y_knots = spline_struct.y

z_knots = evaluate_bicubicspline_interpolant(spline_func, x_knots, y_knots)

plot_topography_with_interpolation_knots(x_int_pts, y_int_pts, z_int_pts,
                                         x_knots, y_knots, z_knots;
                                         xlabel = "x\n [m]",
                                         ylabel = "y\n [m]",
                                         zlabel = "z\n [m]",
                                         azimuth_angle = 54 * pi / 180,
                                         elevation_angle = 27 * pi / 180)

equations = ShallowWaterEquations2D(gravity_constant = 9.81, H0 = 65.0)

function initial_condition_wave(x, t, equations::ShallowWaterEquations2D)
    inicenter = SVector(0.0, 0.0)
    x_norm = x - inicenter
    r = sqrt(x_norm[1]^2 + x_norm[2]^2)

    # Calculate primitive variables
    H = r < 50 ? 75.0 : 65.0
    v1 = 0.0
    v2 = 0.0

    x1, x2 = x
    b = spline_func(x1, x2)

    return prim2cons(SVector(H, v1, v2, b), equations)
end

# Setting initial condition
initial_condition = initial_condition_wave

# Setting the boundary to be a free-slip wall
boundary_condition = boundary_condition_slip_wall

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
solver = DGSEM(polydeg = 3,
               surface_flux = (flux_fjordholm_etal, flux_nonconservative_fjordholm_etal),
               volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

###############################################################################
# Get the TreeMesh and setup a periodic mesh

coordinates_min = (spline_struct.x[1], spline_struct.y[1])
coordinates_max = (spline_struct.x[end], spline_struct.y[end])
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 3,
                n_cells_max = 10_000,
                periodicity = false)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# run the simulation

# define equidistant nodes in time for visualization of an animation
visnodes = range(tspan[1], tspan[2], length = 175)

# use a Runge-Kutta method with error based time step size control
sol = solve(ode, RDPK3SpFSAL49(), abstol = 1.0e-8, reltol = 1.0e-8,
            saveat = visnodes);

# Create an animation of the solution
j = Observable(1)
time = Observable(0.0)

pd_list = [PlotData2D(sol.u[i], semi) for i in 1:length(sol.t)]
f = Figure()

title_string = lift(t -> "time t = $(round(t, digits=3))", time)
az = 130 * pi / 180
el = 18 * pi / 180
ax = Axis3(f[1, 1], xlabel = "x[m]", ylabel = "y[m]", zlabel = "z[m]",
           title = title_string, azimuth = az, elevation = el)

height = lift(i -> pd_list[i].data[1], j)
bottom = lift(i -> pd_list[i].data[4], j)
surface!(ax, pd_list[1].x, pd_list[1].y, bottom;
         colormap = :greenbrownterrain)
wireframe!(ax, pd_list[1].x, pd_list[1].y, height;
           color = Makie.RGBA(0, 0.5, 1, 0.4))
zlims!(ax, 35, 70)

record(f, "animation_2d.gif", 1:length(pd_list)) do tt
    j[] = tt
    time[] = sol.t[tt]
end
