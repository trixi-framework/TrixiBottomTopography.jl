##############################################################
# This script creates data from GeophysicalModelGenerator    #
# into data files which can be used by TrixiBottomTopography #
##############################################################
using GeophysicalModelGenerator
using TrixiBottomTopography
using GMT # needs to be installed for the functions geo_topo_impression and create_topography_data
using DataFrames

# Get a first impression of the topography data
Topo, p, Topo_Cart = geo_topo_impression(resolution = "@earth_relief_01s",
                                         lon_min = 6.963880,
                                         lon_max = 6.978499,
                                         lat_min = 50.947861,
                                         lat_max = 50.957095)

#Create a file with the topography data
df_xyz, Topo_Cart_orth = create_topography_data(low_x = -0.5,
                                                high_x = 0.499,
                                                gridsize_x = 0.001,
                                                low_y = -0.5,
                                                high_y = 0.499,
                                                gridsize_y = 0.002,
                                                write_path = joinpath(@__DIR__, "data"),
                                                dataname = "geo.xyz",
                                                Topo = Topo,
                                                p = p)

# Define file paths
data_dir = joinpath(@__DIR__, "data")
path_src_file = joinpath(data_dir, "geo.xyz")

path_out_file_1d_x = joinpath(data_dir, "rhine_data_1d_20_x_geo.txt")
path_out_file_1d_y = joinpath(data_dir, "rhine_data_1d_20_y_geo.txt")
path_out_file_2d = joinpath(data_dir, "rhine_data_2d_20_geo.txt")

# Convert data
nx = size(Topo_Cart_orth.x.val, 1)
ny = size(Topo_Cart_orth.y.val, 2)
convert_geo_1d(path_src_file,
               path_out_file_1d_x,
               nx = nx,
               ny = ny;
               excerpt = 20,
               section = 10,)
convert_geo_1d(path_src_file,
               path_out_file_1d_y,
               nx = nx,
               ny = ny;
               excerpt = 20,
               direction = "y",
               section = 100,)
convert_geo_2d(path_src_file, path_out_file_2d, nx = nx, ny = ny; excerpt = 20)