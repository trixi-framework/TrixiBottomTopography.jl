##############################################################################
# This script downloads real topography data of the Cliffs of Moher with     #
# GeophysicalModelGenerator.jl and converts it into data files which can be  #
# used by TrixiBottomTopography                                              #
##############################################################################

# Include packages
using TrixiBottomTopography
using GeophysicalModelGenerator
using GMT # needed for the functions `geo_topo_impression` and `create_topography_data`
using DataFrames # needed for the function `create_topography_data`

# Define file paths
root_dir = pkgdir(TrixiBottomTopography)
data_dir = joinpath(root_dir, "examples", "data")

# Get a first impression of the topography data around the Cliffs of Moher, Ireland
Topo, p, Topo_Cart = geo_topo_impression(resolution = "@earth_relief_01s",
                                         lon_min = -9.441139,
                                         lon_max = -9.425224,
                                         lat_min = 52.972167,
                                         lat_max = 52.980888)

# Project the topography onto a regular Cartesian grid and write it to an `xyz` file.
# The bounds `low_x`, `high_x`, `low_y` and `high_y` are given in kilometers and are
# chosen based on the output of `Topo_Cart`.
df_xyz, Topo_Cart_orth = create_topography_data(low_x = -0.55,
                                                high_x = 0.55,
                                                gridsize_x = 0.001,
                                                low_y = -0.55,
                                                high_y = 0.55,
                                                gridsize_y = 0.001,
                                                write_path = data_dir,
                                                dataname = "cliffs.xyz",
                                                Topo = Topo,
                                                p = p)

path_src_file = joinpath(data_dir, "cliffs.xyz")
path_out_file_1d_x = joinpath(data_dir, "cliffs_data_1d_10_x.txt")
path_out_file_2d = joinpath(data_dir, "cliffs_data_2d_10.txt")

# In contrast to `convert_dgm_1d` and `convert_dgm_2d`, the grid dimensions cannot be
# deduced from the `xyz` file itself and therefore have to be passed explicitly
nx = size(Topo_Cart_orth.x.val, 1)
ny = size(Topo_Cart_orth.y.val, 2)

# Convert data. The one dimensional cut is taken in `x` direction through the middle of
# the domain, which runs from the open sea onto the cliffs.
convert_geo_1d(path_src_file, path_out_file_1d_x; nx = nx, ny = ny, excerpt = 10,
               section = 551)
convert_geo_2d(path_src_file, path_out_file_2d; nx = nx, ny = ny, excerpt = 10)
