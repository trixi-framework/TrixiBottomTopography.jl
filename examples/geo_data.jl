#the goal is to use the GeophysicalModelGenerator.jl to load real topography data and use it for TrixiBottomTopography.jl
# 
# the data from originaly used in TrixiBottomTopography.jl  were collected by the Geobasis NRW.
# in this data the first column provides the corresponding ETRS89 East coordinates, the second 
# column the ETRS89 North coordinates and the third column the DHHN2016 height.
# 
# this data is not available in the GeophysicalModelGenerator.jl package, so here we chose first longitudianal
#and transversal data and then choose a projection point in the middle of the choosen area and then convert the data into cartesion 
#data
#
import Pkg
Pkg.add([
            "GeophysicalModelGenerator",
            "GMT",
            "Plots",
            "CSV",
            "DataFrames",
            "TrixiBottomTopography",
            "Downloads",
            "CairoMakie"
        ])
using GMT
using Plots
using GeophysicalModelGenerator
using CSV
using DataFrames
using TrixiBottomTopography
using Downloads: download
using CairoMakie

##########################
# some topography data from the Rhine close to the Theodor-Heuss-Brücke in Mainz
#50.010910 8.275880
#50.009480 8.270871
#50.006184 8.274388
#50.007914 8.278718
###########################
# Specify the limits of the topography data based on the coordinates
lon_min = 8.270871
lon_max = 8.278718
lat_min = 50.006184
lat_max = 50.010910
limits = [lon_min, lon_max, lat_min, lat_max]

lon_mean = (lon_max + lon_min) / 2
lat_mean = (lat_min + lat_max) / 2

###################################
#Loading the topography data
#important: if you have a small area you might use hiher resolution: 
# Note: 
# ====
# - latitude values in the southern hemisphere should have a minus sign (e.g., -2.8)
# - longitude values that are "west" should *either* come with a minus sign *or* are defined by values >180

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

###################################

###
#Lon=17.3, Lat=37.5
#Topo = import_topo(limits, file="@earth_relief_20m")
Topo = import_topo(lon = [8.270871, 8.278718], lat = [50.006184, 50.010910],
                   file = "@earth_relief_01m")
p = ProjectionPoint(Lon = lon_mean, Lat = lat_mean)
Topo_Cart = convert2CartData(Topo, p) # here we get a first impression on what intervall to chose

# the gridpoints have to fullfill the condition: Int(sqrt(length))

low_x = -0.5
high_x = 0.5
gridsize_x = 0.1

low_y = -0.5
high_y = 0.5
gridzize_y = 0.1

values_x = collect(low_x:gridsize_x:high_x)
values_y = collect(low_y:gridzize_y:high_y)

##############
#check if we have a right choice for our gridpoints
function safe_computation(values_x, values_y)
    try
        Int(sqrt(length(values_x) * length(values_y)))
        return Int(sqrt(length(values_x) * length(values_y))), false
    catch e
        if e isa InexactError
            println("There is an InexactError. the gridsizez has to be adjusted")
            return nothing, nothing, true
        else
            rethrow(e)
        end
    end
end
##################

safe_computation(values_x, values_y) # here we check if the gridpoints are ok

Topo_Cart_orth = CartData(xyz_grid(low_x:gridsize_x:high_x, low_y:gridzize_y:high_y, 0))

Topo_Cart_orth = CartData(xyz_grid(-0.5:0.1:0.5, -0.5:0.1:0.5, 0))

Topo_Cart_orth = project_CartData(Topo_Cart_orth, Topo, p)

df_x = DataFrame(Topo_Cart_orth.x.val[:, :, 1], :auto)

df_y = DataFrame(Topo_Cart_orth.y.val[:, :, 1], :auto);

df_z = DataFrame(Topo_Cart_orth.z.val[:, :, 1], :auto);

# Kombiniere die DataFrames für x, y und z in einen DataFrame
df_xyz = DataFrame(x = convert.(Float64, vec(Topo_Cart_orth.x.val[:, :, 1])),  # here we have to convert the values to Float64
                   y = convert.(Float64, vec(Topo_Cart_orth.y.val[:, :, 1])),
                   z = convert.(Float64, vec(Topo_Cart_orth.z.val[:, :, 1])))

# write the data without header and space as delimiter 
# Make sure the data directory exists
data_dir = joinpath(@__DIR__, "data")
mkpath(data_dir)  # Create the directory if it doesn't exist

# Write directly to the data directory
output_file = joinpath(data_dir, "test.xyz")
open(output_file, "w") do file
    for row in eachrow(df_xyz)
        # Round each value in the row to 2 decimal places
        rounded_row = [round(value, digits = 4) for value in row]
        println(file, join(rounded_row, " "))
    end
end

# Download the raw bottom topography data
path_src_file = joinpath(@__DIR__, "test.xyz")

path_out_file_1d_x = joinpath(data_dir, "rhine_data_1d_x_theodor.txt")
path_out_file_1d_y = joinpath(data_dir, "rhine_data_1d_20_y_theodor.txt")
path_out_file_2d = joinpath(data_dir, "rhine_data_2d_20_theodor.txt")

# Convert data
convert_dgm_1d(path_src_file, path_out_file_1d_x; excerpt = 1, section = 1)
convert_dgm_1d(path_src_file, path_out_file_1d_y; excerpt = 1, direction = "y", section = 1)
convert_dgm_2d(path_src_file, path_out_file_2d; excerpt = 1)

#################
#now redoo the steps as in the other turoial steps. try if we can do b-spline interpolation and simulations with this data

# Define data path

data = joinpath(data_dir, "rhine_data_1d_x_theodor.txt")

# Define B-spline structure
spline_struct = CubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 999)

spline_func(x) = spline_interpolation(spline_struct, x)

# Define interpolation points
n = 200
x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))

# Get interpolated values
y_int_pts = spline_func.(x_int_pts)

plot_topography(x_int_pts, y_int_pts; xlabel = "ETRS89 East", ylabel = "DHHN2016 Height")

# Get the original interpolation knots
x_knots = spline_struct.x
y_knots = spline_func.(x_knots)

plot_topography_with_interpolation_knots(x_int_pts, y_int_pts, x_knots, y_knots;
                                         xlabel = "ETRS89 East", ylabel = "DHHN2016 Height")

data = joinpath(data_dir, "rhine_data_2d_20_theodor.txt")

spline_struct = BicubicBSpline(data; end_condition = "not-a-knot", smoothing_factor = 9999)

# Define B-spline interpolation function
spline_func(x, y) = spline_interpolation(spline_struct, x, y)

# Define interpolation points
n = 100
x_int_pts = Vector(LinRange(spline_struct.x[1], spline_struct.x[end], n))
y_int_pts = Vector(LinRange(spline_struct.y[1], spline_struct.y[end], n))

# Get interpolated matrix
z_int_pts = evaluate_bicubicspline_interpolant(spline_func, x_int_pts, y_int_pts)

plot_topography(x_int_pts, y_int_pts, z_int_pts;
                xlabel = "ETRS89\n East",
                ylabel = "ETRS89\n North",
                zlabel = "DHHN2016\n Height",
                azimuth_angle = 54 * pi / 180,
                elevation_angle = 27 * pi / 180)

# Get the original interpolation knots
x_knots = spline_struct.x
y_knots = spline_struct.y
z_knots = evaluate_bicubicspline_interpolant(spline_func, x_knots, y_knots)

plot_topography_with_interpolation_knots(x_int_pts, y_int_pts, z_int_pts,
                                         x_knots, y_knots, z_knots;
                                         xlabel = "ETRS89\n East",
                                         ylabel = "ETRS89\n North",
                                         zlabel = "DHHN2016\n Height",
                                         azimuth_angle = 54 * pi / 180,
                                         elevation_angle = 27 * pi / 180)